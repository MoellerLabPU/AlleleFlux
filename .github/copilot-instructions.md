# AI Rules for {{AlleleFlux}}

- Always start with the simplest clear solution.
- Suggest optimized or advanced solutions only when justified by significant improvements, clearly explaining trade-offs involved.
- Before writing new code, carefully inspect the codebase to ensure similar logic doesn't already exist.
- When introducing new patterns or refactoring, fully remove outdated or duplicate implementations immediately afterward.
- You are careful to only make changes that are requested or you are confident are well understood and related to the change being requested  
- Avoid making major changes to the patterns and architecture of how a feature works, after it has shown to work well, unless explicitly structured
- Introduce new architectures or technologies only after thoroughly exhausting existing solutions, explicitly justifying why the current approach is insufficient.
- Consistently adhere to clear code structure, formatting, and naming conventions. Follow language-specific style guides (e.g., PEP 8 for Python).
- Write explanatory comments when logic is non-obvious or complex, focusing primarily on why the code was written that way. Avoid unnecessary verbosity or commenting on self-explanatory code.
- Carefully evaluate how your changes affect other modules, methods, or workflows.
- Anticipate and proactively handle edge cases, potential race conditions, and security vulnerabilities whenever applicable.
- Include full code examples with necessary import statements, dependencies, and initialization code unless instructed otherwise.
- Employ defensive programming techniques, including validation of inputs, clear error handling, and explicit type checking where beneficial.
- Briefly provide the rationale behind your chosen solutions and reference relevant documentation or learning resources where helpful.
- When addressing issues or bugs, explicitly describe the root cause and clearly explain how your solution resolves the problem. Confirm understanding before proceeding with fixes.
- Favor elegant, maintainable, and idiomatic code.
- Assume readability through descriptive naming and clear structure; reserve comments for deeper insights and rationale.