import calcengine
# Quick debug: check if simplifier handles x/(2*x)
r = calcengine.simplify("x/(2*x)")
print(f"simplify(x/(2*x)) = {r}")
r2 = calcengine.simplify("x*cos(x^2)")
print(f"simplify(x*cos(x^2)) = {r2}")
r3 = calcengine.integrate("x*cos(x^2)", "x")
print(f"integrate(x*cos(x^2), x) = {r3}")
