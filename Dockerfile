# Stage 1: Build WASM using pre-installed Rust toolchain
FROM --platform=$BUILDPLATFORM rust:slim-bookworm AS wasm-builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    && rm -rf /var/lib/apt/lists/*

RUN rustup target add wasm32-unknown-unknown
RUN curl --proto '=https' --tlsv1.2 -sSf https://rustwasm.github.io/wasm-pack/installer/init.sh | sh

WORKDIR /app
COPY scattering-core/ scattering-core/
RUN cd scattering-core && wasm-pack build --target web

# Stage 2: Build frontend
FROM --platform=$BUILDPLATFORM node:20-bookworm-slim AS frontend-builder

WORKDIR /app
COPY package.json package-lock.json ./
RUN npm ci --ignore-scripts

COPY . .
COPY --from=wasm-builder /app/scattering-core/pkg/ ./scattering-core/pkg/
RUN npm run build

# Stage 3: Serve
FROM nginx:alpine
COPY --from=frontend-builder /app/dist /usr/share/nginx/html
EXPOSE 80
CMD ["nginx", "-g", "daemon off;"]
