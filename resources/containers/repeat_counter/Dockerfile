FROM rust:1.83-slim-bookworm AS builder

WORKDIR /usr/src/app
COPY . .
RUN cargo install --path .

# --- Final Stage ---
FROM debian:bookworm-slim

# Install pigz for parallel decompression
RUN apt-get update && apt-get install -y pigz && rm -rf /var/lib/apt/lists/*

COPY --from=builder /usr/local/cargo/bin/repeat_counter /usr/local/bin/repeat_counter
RUN chmod +x /usr/local/bin/repeat_counter

CMD ["repeat_counter", "--help"]

#docker build -t repeatquant:1.0.0 . && apptainer build repeatquant.sif docker-daemon://repeatquant:1.0.0
