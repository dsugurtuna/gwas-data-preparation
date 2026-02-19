FROM python:3.12-slim
RUN apt-get update && apt-get install -y --no-install-recommends plink1.9 bcftools && rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY . .
RUN pip install --no-cache-dir -e ".[dev]"
CMD ["pytest", "-v"]
