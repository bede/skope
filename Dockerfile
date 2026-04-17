FROM rust:1-slim-bookworm AS builder

WORKDIR /build

RUN apt-get update && apt-get install -y --no-install-recommends pkg-config \
    && rm -rf /var/lib/apt/lists/*

ARG TARGETARCH

RUN case "$TARGETARCH" in \
        amd64) printf 'export CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_RUSTFLAGS="-C target-cpu=x86-64-v2"\nexport CARGO_FEATURES="--features scalar"\n' ;; \
        arm64) printf 'export CARGO_TARGET_AARCH64_UNKNOWN_LINUX_GNU_RUSTFLAGS="-C target-cpu=neoverse-n1"\nexport CARGO_FEATURES=""\n' ;; \
        *)     printf 'export CARGO_FEATURES=""\n' ;; \
    esac > /cargo_rustflags.sh

COPY Cargo.toml Cargo.lock ./
RUN mkdir src && echo "pub fn main() {}" > src/main.rs && echo > src/lib.rs
RUN . /cargo_rustflags.sh && cargo build --release --locked $CARGO_FEATURES 2>/dev/null; true

COPY src ./src
RUN . /cargo_rustflags.sh && touch src/main.rs src/lib.rs && cargo build --release --locked $CARGO_FEATURES

FROM debian:bookworm-slim

RUN apt-get update && apt-get install -y --no-install-recommends bash \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /build/target/release/skope /usr/local/bin/skope

ENTRYPOINT ["/bin/bash", "-c"]
CMD ["skope"]
