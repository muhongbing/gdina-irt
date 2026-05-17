# plumber.R - 精准诊断 API（无视觉/大模型版）
library(plumber)
library(mirt)
library(GDINA)
library(lme4)
library(jsonlite)
library(dplyr)
library(tidyr)

# 尝试加载 gtheory，失败则回退到 lme4
use_gtheory <- suppressWarnings(suppressMessages({
  if (!requireNamespace("gtheory", quietly = TRUE)) {
    tryCatch({
      remotes::install_version("gtheory", version = "0.1.2", quiet = TRUE)
      library(gtheory, quietly = TRUE)
      TRUE
    }, error = function(e) FALSE)
  } else {
    library(gtheory, quietly = TRUE)
    TRUE
  }
}))

# ======================== 辅助函数 ========================

arbitrate_scores <- function(scores_a, scores_b, n_items) {
  if (is.null(scores_a) && is.null(scores_b)) {
    return(list(final_scores = rep(0, n_items),
                conflict_count = 0,
                flags = rep("both_failed", n_items),
                error = TRUE))
  }
  if (is.null(scores_a)) {
    final <- sapply(scores_b, function(x) x$score)
    return(list(final_scores = final, conflict_count = 0,
                flags = rep("single", length(final)), error = FALSE))
  }
  if (is.null(scores_b)) {
    final <- sapply(scores_a, function(x) x$score)
    return(list(final_scores = final, conflict_count = 0,
                flags = rep("single", length(final)), error = FALSE))
  }
  n <- max(length(scores_a), length(scores_b))
  final <- numeric(n)
  flags <- character(n)
  for (i in seq_len(n)) {
    sa <- if (i <= length(scores_a)) scores_a[[i]]$score else NA
    sb <- if (i <= length(scores_b)) scores_b[[i]]$score else NA
    if (is.na(sa) && is.na(sb)) {
      final[i] <- 0
      flags[i] <- "both_missing"
    } else if (is.na(sa)) {
      final[i] <- sb
      flags[i] <- "single"
    } else if (is.na(sb)) {
      final[i] <- sa
      flags[i] <- "single"
    } else {
      diff <- abs(sa - sb)
      if (diff == 0) {
        final[i] <- sa
        flags[i] <- "agree"
      } else if (diff <= 1) {
        final[i] <- round((sa + sb) / 2)
        flags[i] <- "minor_diff"
      } else {
        final[i] <- sa
        flags[i] <- "conflict"
      }
    }
  }
  list(final_scores = final, conflict_count = sum(flags == "conflict"),
       flags = flags, error = FALSE)
}

discretize_scores <- function(raw_scores, rubric_df) {
  discrete <- numeric(length(raw_scores))
  for (i in seq_along(raw_scores)) {
    score <- raw_scores[i]
    if (i > nrow(rubric_df)) {
      discrete[i] <- score
      next
    }
    thr <- rubric_df$thresholds[[i]]
    max_s <- rubric_df$max_scores[i]
    if (!is.null(thr) && length(thr) > 0) {
      discrete[i] <- sum(score >= thr)
    } else if (max_s <= 1) {
      discrete[i] <- round(score)
    } else {
      discrete[i] <- pmax(0, round(score))
    }
  }
  discrete
}

run_diagnosis <- function(student_ids, scores_m1, scores_m2,
                          rubric_df, Q_list,
                          weights = c(0.4, 0.6),
                          sequential = TRUE) {
  stopifnot(is.matrix(scores_m1), is.matrix(scores_m2))
  stopifnot(nrow(scores_m1) == nrow(scores_m2),
            ncol(scores_m1) == ncol(scores_m2))
  N <- nrow(scores_m1)
  J <- ncol(scores_m1)
  stopifnot(length(student_ids) == N)
  stopifnot(nrow(rubric_df) == J)

  # 仲裁两个评分员
  score_list_a <- lapply(seq_len(N), function(i) {
    lapply(seq_len(J), function(j) list(score = scores_m1[i, j]))
  })
  score_list_b <- lapply(seq_len(N), function(i) {
    lapply(seq_len(J), function(j) list(score = scores_m2[i, j]))
  })
  arb_results <- lapply(seq_len(N), function(i) {
    arbitrate_scores(score_list_a[[i]], score_list_b[[i]], J)
  })
  combined_raw <- do.call(rbind, lapply(arb_results, `[[`, "final_scores"))
  flags_mat <- do.call(rbind, lapply(arb_results, `[[`, "flags"))
  colnames(combined_raw) <- rubric_df$item_id
  rownames(combined_raw) <- student_ids
  colnames(flags_mat) <- paste0("Item", seq_len(J))
  rownames(flags_mat) <- student_ids

  # IRT 估计权重 lambda 及能力 theta
  theta <- rep(0, N)
  lambda <- 0.5
  se_lambda <- NA
  combined <- combined_raw
  tryCatch({
    init <- combined
    zero_var <- which(apply(init, 2, function(x) length(unique(x)) < 2) |
                        apply(init, 2, sd, na.rm = TRUE) < 1e-10)
    if (length(zero_var) > 0) {
      set.seed(42 + N)
      init[, zero_var] <- init[, zero_var] +
        matrix(rnorm(N * length(zero_var), 0, 0.05), N, length(zero_var))
    }
    gpcm <- mirt(init, 1, itemtype = "gpcm", verbose = FALSE, SE = TRUE)
    theta <- fscores(gpcm, full.scores.SE = TRUE)[, 1]
    items <- coef(gpcm, IRTpars = TRUE, simplify = TRUE)$items
    compute_exp <- function(th) {
      sapply(seq_len(J), function(j) {
        a <- items[j, "a"]
        b <- as.numeric(items[j, grep("^b[0-9]+$", colnames(items))])
        b <- b[!is.na(b)]
        sum(plogis(a * (th - b)))
      })
    }
    exp_scores <- t(sapply(theta, compute_exp))
    D <- scores_m1 - scores_m2
    num <- sum(D * (exp_scores - scores_m2), na.rm = TRUE)
    den <- sum(D^2, na.rm = TRUE)
    lambda <- if (den < 1e-12) 0.5 else pmax(0, pmin(1, num / den))
    combined <- lambda * scores_m1 + (1 - lambda) * scores_m2
    se_lambda <- sqrt(mean((combined - exp_scores)^2, na.rm = TRUE) / den)
  }, error = function(e) {
    message("IRT 失败，使用总分 Z 分数: ", e$message)
    total_scores <- rowSums(combined, na.rm = TRUE)
    theta <- scale(total_scores)[, 1]
    theta[is.na(theta)] <- 0
  })

  # 离散化得分
  discrete <- matrix(0, N, J)
  types <- rubric_df$item_type
  thresholds <- rubric_df$thresholds
  for (j in seq_len(J)) {
    if (types[j] == "objective") {
      discrete[, j] <- ifelse(combined[, j] >= 0.5, 1, 0)
    } else {
      thr <- thresholds[[j]]
      score_j <- combined[, j]
      if (is.null(thr)) {
        discrete[, j] <- pmax(0, round(score_j))
      } else {
        breaks <- c(-Inf, thr, Inf)
        level <- as.integer(cut(score_j, breaks = breaks, right = FALSE)) - 1
        level[score_j >= max(thr)] <- length(thr)
        discrete[, j] <- level
      }
    }
  }

  # 多元概化理论信度
  phi <- phi_obj <- phi_subj <- NA
  tryCatch({
    df <- data.frame(
      Person = factor(rep(seq_len(N), each = J * 2)),
      Rater  = factor(rep(rep(c("M1", "M2"), each = J), N)),
      Item   = factor(rep(seq_len(J), N * 2)),
      Type   = factor(rep(ifelse(types == "objective", "obj", "subj"), N * 2)),
      Score  = c(t(scores_m1), t(scores_m2))
    )
    if (use_gtheory) {
      g_obj <- gstudy(subset(df, Type == "obj"),
                      formula = Score ~ (1 | Person) + (1 | Rater) + (1 | Item) +
                        (1 | Person:Rater) + (1 | Person:Item) + (1 | Rater:Item))
      g_subj <- gstudy(subset(df, Type == "subj"),
                       formula = Score ~ (1 | Person) + (1 | Rater) + (1 | Item) +
                         (1 | Person:Rater) + (1 | Person:Item) + (1 | Rater:Item))
      extract_vc <- function(g) {
        comps <- g$components
        get_var <- function(name) if (!is.null(comps[[name]])) comps[[name]]$var else 0
        list(p = get_var("Person"), r = get_var("Rater"), i = get_var("Item"),
             pr = get_var("Person:Rater"), pi = get_var("Person:Item"),
             ri = get_var("Rater:Item"), e = get_var("Residual"))
      }
      vc_obj <- extract_vc(g_obj)
      vc_subj <- extract_vc(g_subj)
    } else {
      fit_obj <- lmer(Score ~ (1 | Person) + (1 | Rater) + (1 | Item) +
                        (1 | Person:Rater) + (1 | Person:Item) + (1 | Rater:Item),
                      data = subset(df, Type == "obj"), REML = TRUE)
      fit_subj <- lmer(Score ~ (1 | Person) + (1 | Rater) + (1 | Item) +
                         (1 | Person:Rater) + (1 | Person:Item) + (1 | Rater:Item),
                       data = subset(df, Type == "subj"), REML = TRUE)
      extract_vc <- function(m) {
        vc <- as.data.frame(VarCorr(m))
        grp <- vc$grp
        get_var <- function(g) {
          idx <- which(grp == g)
          if (length(idx) > 0) max(vc$vcov[idx[1]], 0) else 0
        }
        resid_var <- if ("Residual" %in% grp) get_var("Residual") else attr(VarCorr(m), "sc")^2
        list(p = get_var("Person"), r = get_var("Rater"), i = get_var("Item"),
             pr = get_var("Person:Rater"), pi = get_var("Person:Item"),
             ri = get_var("Rater:Item"), e = resid_var)
      }
      vc_obj <- extract_vc(fit_obj)
      vc_subj <- extract_vc(fit_subj)
    }
    n_obj <- sum(types == "objective")
    n_subj <- J - n_obj
    n_r <- 2
    err_obj <- with(vc_obj, r/n_r + i/n_obj + pr/n_r + pi/n_obj + ri/(n_r * n_obj) + e/(n_r * n_obj))
    err_subj <- with(vc_subj, r/n_r + i/n_subj + pr/n_r + pi/n_subj + ri/(n_r * n_subj) + e/(n_r * n_subj))
    w <- weights
    var_true <- w[1]^2 * vc_obj$p + w[2]^2 * vc_subj$p
    var_error <- w[1]^2 * err_obj + w[2]^2 * err_subj
    phi <- var_true / (var_true + var_error)
    phi_obj <- vc_obj$p / (vc_obj$p + err_obj)
    phi_subj <- vc_subj$p / (vc_subj$p + err_subj)
  }, error = function(e) {
    message("MGT 失败: ", e$message)
  })

  # GDINA 认知诊断
  mastery <- NULL
  suggestions <- NULL
  fit_metrics <- NULL
  if (N >= 5) {
    tryCatch({
      gdina_fit <- GDINA(dat = discrete, Q = Q_list, model = "GDINA",
                         type = "ordinal", sequential = sequential,
                         control = list(maxitr = 1000, verbose = FALSE))
      mastery <- personparm(gdina_fit, what = "mp")
      suggestions <- apply(mastery, 1, function(m) {
        weak <- order(m)[1:min(3, ncol(mastery))]
        paste("加强:", paste(colnames(Q_list[[1]])[weak], collapse = ", "))
      })
      fit_stats <- tryCatch(modelfit(gdina_fit), error = function(e) NULL)
      if (!is.null(fit_stats)) {
        fit_metrics <- list(AIC = fit_stats$AIC, BIC = fit_stats$BIC,
                            SRMSR = fit_stats$SRMSR, MAD = fit_stats$MAD)
      }
    }, error = function(e) {
      message("GDINA 失败: ", e$message)
    })
  }

  list(
    student_ids = student_ids,
    combined = combined,
    discrete = discrete,
    theta = theta,
    lambda = lambda,
    se_lambda = se_lambda,
    phi = phi,
    phi_obj = phi_obj,
    phi_subj = phi_subj,
    mastery = mastery,
    suggestions = suggestions,
    gdina_fit_metrics = fit_metrics,
    conflict_flags = flags_mat,
    n_students = N,
    n_items = J
  )
}

# ======================== API 端点 ========================

#* @apiTitle 精准诊断 API
#* @apiDescription 基于两个评分员得分矩阵、评分标准和 Q 矩阵进行诊断分析
#* @post /diagnose
function(req, res) {
  tryCatch({
    body <- jsonlite::fromJSON(req$postBody, simplifyVector = TRUE, simplifyMatrix = TRUE)
    
    student_ids <- body$student_ids
    scores_m1 <- as.matrix(body$scores_m1)
    scores_m2 <- as.matrix(body$scores_m2)
    rubric_raw <- body$rubric
    Q_list_raw <- body$Q_list
    weights <- if (!is.null(body$weights)) body$weights else c(0.4, 0.6)
    sequential <- if (!is.null(body$sequential)) body$sequential else TRUE

    rubric_df <- data.frame(
      item_id = rubric_raw$item_ids,
      item_type = rubric_raw$item_types,
      max_score = rubric_raw$max_scores,
      thresholds = I(rubric_raw$thresholds),
      correct_answer = rubric_raw$correct_answers,
      stringsAsFactors = FALSE
    )
    rubric_df$thresholds <- lapply(rubric_df$thresholds, function(x) {
      if (is.null(x) || (is.list(x) && length(x) == 0)) NULL else as.numeric(unlist(x))
    })

    Q_list <- lapply(Q_list_raw, as.matrix)
    names(Q_list) <- rubric_df$item_id

    result <- run_diagnosis(student_ids, scores_m1, scores_m2,
                            rubric_df, Q_list, weights, sequential)
    
    result$combined <- as.matrix(result$combined)
    result$discrete <- as.matrix(result$discrete)
    if (!is.null(result$mastery)) result$mastery <- as.matrix(result$mastery)
    result$conflict_flags <- as.matrix(result$conflict_flags)

    list(success = TRUE, data = result)
  }, error = function(e) {
    res$status <- 500
    list(success = FALSE, error = e$message)
  })
}

#* @get /health
function() {
  list(status = "OK", timestamp = Sys.time())
}