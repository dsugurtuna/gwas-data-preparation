FROM python:3.12-slim
# PLINK 1.9 from Debian, for running merges and conversions for real.
RUN apt-get update && apt-get install -y --no-install-recommends plink1.9 \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY . .
RUN pip install --no-cache-dir -e ".[dev]"
CMD ["pytest", "-v"]
