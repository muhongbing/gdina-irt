FROM rocker/r-ver:4.5.1

# 1. 更换为清华 Ubuntu 源（稳定且速度快）
RUN sed -i 's|http://archive.ubuntu.com|https://mirrors.tuna.tsinghua.edu.cn/ubuntu|g' /etc/apt/sources.list && \
    sed -i 's|http://security.ubuntu.com|https://mirrors.tuna.tsinghua.edu.cn/ubuntu|g' /etc/apt/sources.list

# 2. 安装系统依赖
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    make \
    && rm -rf /var/lib/apt/lists/*

# 3. 设置 CRAN 镜像（清华源）和超时时间
RUN echo 'options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"), timeout = 1200, download.file.method = "libcurl")' >> /usr/local/lib/R/etc/Rprofile.site

# 4. 安装 plumber（启用依赖，确保完整安装）
RUN R -e "install.packages('plumber', dependencies = TRUE, Ncpus = 2)"

# 5. 安装其他核心包（不自动安装依赖，避免 arrow 等大包，但 plumber 已独立安装，不会重复）
RUN R -e "install.packages(c('mirt', 'GDINA', 'lme4', 'jsonlite', 'dplyr', 'tidyr', 'remotes'), dependencies = FALSE, Ncpus = 2)"

# 6. 安装 gtheory（使用单核、长超时，防止内存/编译崩溃）
RUN R -e "install.packages('gtheory', dependencies = FALSE, Ncpus = 1, timeout = 1800)"

# 7. 复制你的 API 脚本
COPY plumber.R /plumber.R

# 8. 暴露端口
EXPOSE 8000

# 9. 启动服务
ENTRYPOINT ["R", "-e", "pr <- plumber::plumb('/plumber.R'); pr$run(host='0.0.0.0', port=8000)"]