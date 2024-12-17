# Coding style
The project enforces a consistent coding style through [`clang-format`](https://clang.llvm.org/docs/ClangFormat.html).
The git hooks to conform the code to the style before any commit should be installed by typing `pre-commit install` once
every time the repository is cloned (requires a working installation of [pre-commit](https://pre-commit.com/#install)).

Similarly to what is done by `lb-format` for LHCb, you can run all the `pre-commit` hooks before committing the code
via `pre-commit run --all` or `pre-commit run --files <file paths>` to avoid having to stage the files twice.
