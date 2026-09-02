# Chalk

**A calculus engine written from scratch in C++, served as a website.**

[![Live Demo](https://img.shields.io/badge/demo-chalk--z9gs.onrender.com-brightgreen)](https://chalk-z9gs.onrender.com)
![C++17](https://img.shields.io/badge/C%2B%2B-17-blue)
![Flask](https://img.shields.io/badge/backend-Flask-lightgrey)
![pybind11](https://img.shields.io/badge/bridge-pybind11-orange)

**[Try it live →](https://chalk-z9gs.onrender.com)** *(free-tier host, first load after idle takes ~30s to wake up)*

No SymPy. No Wolfram Alpha API. No external math library. Every derivative,
integral, limit, and equation is computed by a hand-written C++ engine —
its own tokenizer, parser, and symbolic manipulation code — that Python
just calls into.

---

## Why

Most "build a calculator" web projects reach for SymPy and call it done.
Chalk does the opposite: the actual calculus — parsing an expression into
an AST, differentiating it symbolically, pattern-matching integration
rules, stepping through L'Hôpital's rule — is implemented from first
principles in C++. Python's only job is to be a thin Flask wrapper around
a compiled `.so` module.

## How it's built

| Layer | Role |
|---|---|
| **C++17** (`cpp/engine.cpp`) | The actual math: expression parser, symbolic differentiation, pattern-based integration, limits (direct substitution → L'Hôpital → numeric fallback), equation solving (exact for linear/quadratic, Newton's method otherwise). |
| **pybind11** (`cpp/bindings.cpp`) | Exposes the C++ functions to Python, compiled into a `calcengine.so` module. |
| **Flask** (`python_app/app.py`) | A thin wrapper that calls into `calcengine.so`. Does no math itself. |
| **Frontend** (`python_app/templates/index.html`) | One HTML page, five tools: Differentiate, Integrate, Limit, Solve, Evaluate. |

```
calc_project/
├── cpp/
│   ├── engine.hpp / engine.cpp   # the actual math (parser, calc, solver)
│   ├── bindings.cpp              # pybind11 glue
│   └── test_main.cpp             # standalone C++ test harness
├── python_app/
│   ├── app.py                    # Flask routes -> calcengine
│   └── templates/index.html      # the web UI
├── build.sh                      # compiles cpp/ into python_app/calcengine.so
├── Dockerfile                    # for deploying as a container
└── requirements.txt
```

## What the engine can and can't do

**Differentiation** — full support: products, quotients, chain rule, all
standard functions (`sin cos tan asin acos atan ln log exp sqrt abs`), and
the generalized power rule (handles `x^x`-style expressions too).

- Evaluate symbolic vectors and paths
- **Ordinary Differential Equations (ODE)**
  - Solve 1st-order and simple 2nd-order ODEs
  - Interactive slope fields and numerical solution trajectories (RK4)
  > **Note**: The ODE solver is currently slow and occasionally unstable (crashing). This will be fixed in a later update.
- **Series & Sequences**
  - Taylor series expansion up to a specified order.

**Integration** — polynomials, `1/x`, constant multiples/sums, and
linear-argument trig/exp/ln patterns (e.g. `sin(3x+1)`, `1/(2x+3)`,
`exp(2x)`). No integration by parts or trig substitution — expressions
must match known atomic forms.
- **New Feature**: Supports indefinite multiple integrals (e.g. `x, y` for multivariable integration).

**Limits** — evaluates $\lim_{x \to c} f(x)$ using direct substitution and L'Hôpital's rule
or $\infty / \infty$ forms (recursively, up to depth 6), then a numeric two-sided
fallback. Handles limits at `±infinity` too.

**Special Functions** — native parser and evaluation support for Error function (`erf(x)`) and Fresnel Integrals (`S(x)`, `C(x)` or `FresnelS(x)`, `FresnelC(x)`).

**Syntax Guide** — click the `i` info button in the top left of the UI for an interactive syntax guide detailing all supported math symbols, functions, and constants.

**Equations** — exact solutions for linear and quadratic equations
(including complex roots), detected automatically by checking whether the
third derivative vanishes. Higher-degree equations fall back to Newton's
method seeded across `[-20, 20]`.

## Run it locally

```bash
./build.sh                 # compiles the C++ engine into python_app/
cd python_app
python3 app.py             # http://localhost:5000
```



## Roadmap

- [ ] Integration by parts and trig substitution
- [ ] Multi-variable calculus (partial derivatives, gradients)
- [ ] Step-by-step solution output, not just final answers
- [ ] Expression history / shareable permalinks

## License

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.
