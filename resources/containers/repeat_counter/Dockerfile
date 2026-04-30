FROM rust:1.85.0-slim-bookworm AS builder

WORKDIR /usr/src/app
COPY . .
RUN cargo install --locked --path .

# --- Final Stage ---
FROM debian:bookworm-slim

RUN apt-get update && apt-get install -y pigz && rm -rf /var/lib/apt/lists/*

COPY --from=builder /usr/local/cargo/bin/repeat_counter /usr/local/bin/repeat_counter
RUN chmod +x /usr/local/bin/repeat_counter

CMD ["repeat_counter", "--help"]

#docker build -t repeatquant:1.0.1 .
