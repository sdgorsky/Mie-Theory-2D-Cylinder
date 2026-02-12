.PHONY: setup run build build-wasm clean format test docker help

# Colors for output
RED := \033[0;31m
GREEN := \033[0;32m
YELLOW := \033[0;33m
NC := \033[0m # No Color

help:
	@echo "Available commands:"
	@echo "  make setup    - Check and install dependencies"
	@echo "  make run      - Build and run the development server"
	@echo "  make build    - Build WASM and frontend"
	@echo "  make format   - Format and lint Rust and TypeScript code"
	@echo "  make test     - Run tests"
	@echo "  make docker   - Build and run in Docker"
	@echo "  make clean    - Clean build artifacts"

setup:
	@echo "Checking and installing dependencies..."
	@echo ""
	@FAIL=0; \
	if [ -n "$$CI" ]; then AUTO=y; else AUTO=; fi; \
	if ! command -v node &> /dev/null; then \
		echo "$(RED)✗ node is not installed$(NC)"; \
		echo "  Install Node.js from https://nodejs.org/ or via your package manager"; \
		FAIL=1; \
	else \
		echo "$(GREEN)✓ node $$(node --version)$(NC)"; \
	fi; \
	if ! command -v npm &> /dev/null; then \
		echo "$(RED)✗ npm is not installed$(NC)"; \
		echo "  npm is included with Node.js — install Node first"; \
		FAIL=1; \
	else \
		echo "$(GREEN)✓ npm $$(npm --version)$(NC)"; \
	fi; \
	if ! command -v rustup &> /dev/null; then \
		echo "$(YELLOW)✗ Rust is not installed$(NC)"; \
		if [ -n "$$AUTO" ]; then confirm=y; else read -p "  Install Rust via rustup? [y/N] " confirm; fi; \
		if [ "$$confirm" = "y" ] || [ "$$confirm" = "Y" ]; then \
			curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y; \
			. "$$HOME/.cargo/env"; \
			echo "$(GREEN)✓ Rust installed ($$(cargo --version))$(NC)"; \
		else \
			FAIL=1; \
		fi; \
	else \
		echo "$(GREEN)✓ cargo $$(cargo --version 2>/dev/null | cut -d' ' -f2)$(NC)"; \
	fi; \
	if command -v rustup &> /dev/null; then \
		if ! rustup target list --installed 2>/dev/null | grep -q wasm32-unknown-unknown; then \
			echo "$(YELLOW)✗ wasm32-unknown-unknown target missing$(NC)"; \
			if [ -n "$$AUTO" ]; then confirm=y; else read -p "  Install wasm32 target? [y/N] " confirm; fi; \
			if [ "$$confirm" = "y" ] || [ "$$confirm" = "Y" ]; then \
				rustup target add wasm32-unknown-unknown; \
				echo "$(GREEN)✓ wasm32-unknown-unknown target installed$(NC)"; \
			else \
				FAIL=1; \
			fi; \
		else \
			echo "$(GREEN)✓ wasm32-unknown-unknown target$(NC)"; \
		fi; \
		if ! command -v wasm-pack &> /dev/null; then \
			echo "$(YELLOW)✗ wasm-pack is not installed$(NC)"; \
			if [ -n "$$AUTO" ]; then confirm=y; else read -p "  Install wasm-pack? [y/N] " confirm; fi; \
			if [ "$$confirm" = "y" ] || [ "$$confirm" = "Y" ]; then \
				cargo install wasm-pack; \
				echo "$(GREEN)✓ wasm-pack installed$(NC)"; \
			else \
				FAIL=1; \
			fi; \
		else \
			echo "$(GREEN)✓ wasm-pack $$(wasm-pack --version 2>/dev/null)$(NC)"; \
		fi; \
		if ! command -v rustfmt &> /dev/null; then \
			echo "$(YELLOW)✗ rustfmt is not installed$(NC)"; \
			if [ -n "$$AUTO" ]; then confirm=y; else read -p "  Install rustfmt? [y/N] " confirm; fi; \
			if [ "$$confirm" = "y" ] || [ "$$confirm" = "Y" ]; then \
				rustup component add rustfmt; \
				echo "$(GREEN)✓ rustfmt installed$(NC)"; \
			else \
				FAIL=1; \
			fi; \
		else \
			echo "$(GREEN)✓ rustfmt$(NC)"; \
		fi; \
		if ! rustup component list 2>/dev/null | grep -q "clippy.*installed"; then \
			echo "$(YELLOW)✗ clippy is not installed$(NC)"; \
			if [ -n "$$AUTO" ]; then confirm=y; else read -p "  Install clippy? [y/N] " confirm; fi; \
			if [ "$$confirm" = "y" ] || [ "$$confirm" = "Y" ]; then \
				rustup component add clippy; \
				echo "$(GREEN)✓ clippy installed$(NC)"; \
			else \
				FAIL=1; \
			fi; \
		else \
			echo "$(GREEN)✓ clippy$(NC)"; \
		fi; \
	fi; \
	if command -v npm &> /dev/null && [ ! -d "node_modules" ]; then \
		echo "$(YELLOW)✗ node_modules not installed$(NC)"; \
		if [ -n "$$AUTO" ]; then confirm=y; else read -p "  Run npm install? [y/N] " confirm; fi; \
		if [ "$$confirm" = "y" ] || [ "$$confirm" = "Y" ]; then \
			npm install; \
			echo "$(GREEN)✓ npm dependencies installed$(NC)"; \
		else \
			FAIL=1; \
		fi; \
	elif [ -d "node_modules" ]; then \
		echo "$(GREEN)✓ node_modules$(NC)"; \
	fi; \
	echo ""; \
	if [ "$$FAIL" = "1" ]; then \
		echo "$(YELLOW)Some dependencies are missing. Re-run 'make setup' after installing them.$(NC)"; \
		exit 1; \
	fi
	@git config core.hooksPath .githooks
	@echo "$(GREEN)Git hooks configured.$(NC)"
	@echo ""
	@echo "$(GREEN)Setup complete! Run 'make run' to start.$(NC)"

# Build WASM module
build-wasm:
	@echo "Building WASM module..."
	cd scattering-core && wasm-pack build --target web


# Full production build
build: build-wasm
	@if [ ! -d "node_modules" ]; then \
		echo "Installing npm dependencies..."; \
		npm install; \
	fi
	@echo "Building frontend..."
	npm run build
	@echo "$(GREEN)Build complete!$(NC)"

# Run development server
run: build-wasm
	@if [ ! -d "node_modules" ]; then \
		echo "Installing npm dependencies..."; \
		npm install; \
	fi
	@echo "Starting development server..."
	npm run dev
	

# Format and lint both Rust and TypeScript
format:
	@echo "Formatting Rust code..."
	cd scattering-core && cargo fmt
	@echo "Linting Rust code..."
	cd scattering-core && cargo clippy -- -D warnings
	@echo "Formatting TypeScript code..."
	npm run format
	@echo "Linting TypeScript code..."
	npm run lint
	@echo "$(GREEN)Formatting and linting complete!$(NC)"

# Run all tests
test:
	@echo "Running Rust tests..."
	cd scattering-core && cargo test
	@echo "$(GREEN)All tests passed!$(NC)"

# Build and run in Docker
docker:
	@echo "Building Docker image..."
	docker build -t 2d-scattering .
	@echo "$(GREEN)Docker image built!$(NC)"
	@echo "Starting container at http://localhost:8080..."
	docker run --rm -p 8080:80 2d-scattering

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
	cd scattering-core && cargo clean
	rm -rf scattering-core/pkg
	rm -rf node_modules
	rm -rf dist
	@echo "$(GREEN)Clean complete!$(NC)"
