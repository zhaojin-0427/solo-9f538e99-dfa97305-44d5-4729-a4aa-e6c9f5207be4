FROM python:3.12-slim

WORKDIR /srv/app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

ENV APP_DATABASE=/data/app.db \
    COMPUTE_DELAY_SECONDS=2 \
    PORT=5000

RUN mkdir -p /data
VOLUME ["/data"]
EXPOSE 5000

CMD ["python", "run.py"]
