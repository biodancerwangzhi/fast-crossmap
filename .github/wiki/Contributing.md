# Contributing

Thank you for considering contributing to FastCrossMap!

## Ways to Contribute

- Report bugs
- Suggest features
- Improve documentation
- Submit code changes

---

## Reporting Issues

When reporting bugs, please include:

1. FastCrossMap version: `fast-crossmap --version`
2. Operating system and version
3. Command that caused the issue
4. Error message (if any)
5. Sample data (if possible)

---

## Code Contributions

### Development Setup

```bash
# 1. Install Rust (1.70+)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env

# 2. Clone and build
git clone https://github.com/YOUR_USERNAME/fast-crossmap.git
cd fast-crossmap
cargo build --release

# 3. Run tests
cargo test
```

### Development Workflow

```bash
# Create feature branch
git checkout -b feature/my-feature

# Make changes and test
cargo build
cargo test

# Format and lint
cargo fmt
cargo clippy

# Commit and push
git add .
git commit -m "Add my feature"
git push origin feature/my-feature
```

Then open a Pull Request on GitHub.

### Code Style

- Follow Rust conventions
- Run `cargo fmt` before committing
- Run `cargo clippy` to check for issues
- Add tests for new features

### Testing

```bash
cargo test              # Run all tests
cargo test test_name    # Run specific test
cargo test -- --nocapture  # With output
cargo bench             # Run benchmarks
```

### Project Structure

```
fast-crossmap/
├── src/
│   ├── main.rs          # CLI entry point
│   ├── lib.rs           # Library root
│   ├── formats/         # Format implementations
│   │   ├── bed.rs
│   │   ├── bam.rs
│   │   ├── vcf.rs
│   │   └── ...
│   └── ...
├── tests/               # Integration tests
├── benches/             # Benchmarks
└── Cargo.toml           # Dependencies
```

---

## Documentation

Help improve documentation:
- Fix typos
- Add examples
- Clarify explanations

---

## License

By contributing, you agree that your contributions will be licensed under the MIT License.

---

## Contact

- GitHub Issues: https://github.com/biodancerwangzhi/fast-crossmap/issues
- Email: szxszx@foxmail.com
