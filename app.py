

# ============================================================
# CONFIGURAÇÃO DA PÁGINA
# ============================================================
st.set_page_config(page_title="Métodos Numéricos - Raízes", layout="wide")


# ============================================================
# FUNÇÃO PARA PROCESSAR EXPRESSÃO MATEMÁTICA
# ============================================================
def criar_funcao(expressao):
    expr = expressao.strip()
    expr = expr.replace("^", "**")
    expr = expr.replace("X", "x")
    expr = re.sub(r"(\d+),(\d+)", r"\1.\2", expr)

    substituicoes = {
        "sen": "sin",
        "tg": "tan",
        "raiz": "sqrt",
        "ln": "log",
    }
    for a, b in substituicoes.items():
        expr = re.sub(rf"\b{a}\b", b, expr, flags=re.IGNORECASE)

    # Multiplicação implícita 2x, 3(x+1), x(x+1), )(, 2pi, xpi
    expr = re.sub(r"(\d+)([a-zA-Z])", r"\1*\2", expr)
    expr = re.sub(r"(\d+)(\()", r"\1*\2", expr)
    expr = re.sub(r"([a-zA-Z])(\()", r"\1*\2", expr) expr = re.sub(r"(\))(\()", r"\1*\2", expr) expr = re.sub(r"(\))([a-zA-Z])", r"\1*\2", expr)

    def f(x):
        ambiente = {
            "x": x,
            "e": math.e,
            "pi": math.pi,
            "exp": math.exp,
            "sin": math.sin,
            "cos": math.cos,
            "tan": math.tan,
            "sqrt": math.sqrt,
            "log": math.log,
            "abs": abs,
            "pow": pow,
        }
        return float(eval(expr, {"__builtins__": {}}, ambiente))

    # Testa a expressão no ponto 0.0 para garantir validade sintática
    try:
        f(0.0)
    except Exception:
        pass

    return f, expr


# ============================================================
# MÉTODOS NUMÉRICOS
# ============================================================
def metodo_bissecao(f, a, b, tol, max_iter):
    fa, fb = f(a), f(b)
    if fa * fb > 0:
        raise ValueError("Na Bisseção, f(a) e f(b) precisam ter sinais opostos.")

    linhas = []
    geometria = []

    for i in range(1, max_iter + 1):
        c = (a + b) / 2.0
        fc = f(c)
        erro = abs(b - a) / 2.0

        linhas.append([i, a, b, c, fa, fb, fc, erro])
        geometria.append((a, b, c))

        if abs(fc) < tol or erro < tol:
            return (
                c,
                pd.DataFrame(
                    linhas,
                    columns=["Iteração", "a", "b", "c", "f(a)", "f(b)", "f(c)", "Erro"],
                ),
                geometria,
            )

        if fa * fc < 0:
            b, fb = c, fc
        else:
            a, fa = c, fc

    return (
        c,
        pd.DataFrame(
            linhas,
            columns=["Iteração", "a", "b", "c", "f(a)", "f(b)", "f(c)", "Erro"],
        ),
        geometria,
    )


def metodo_falsa_posicao(f, a, b, tol, max_iter):
    fa, fb = f(a), f(b)
    if fa * fb > 0:
        raise ValueError("Na Falsa Posição, f(a) e f(b) precisam ter sinais opostos.")

    linhas = []
    geometria = []

    for i in range(1, max_iter + 1):
        c = (a * fb - b * fa) / (fb - fa)
        fc = f(c)
        erro = abs(fc)

        linhas.append([i, a, b, c, fa, fb, fc, erro])
        geometria.append((a, fa, b, fb, c, fc))

        if abs(fc) < tol:
            return (
                c,
                pd.DataFrame(
                    linhas,
                    columns=["Iteração", "a", "b", "c", "f(a)", "f(b)", "f(c)", "Erro"],
                ),
                geometria,
            )

        if fa * fc < 0:
            b, fb = c, fc
        else:
            a, fa = c, fc

    return (
        c,
        pd.DataFrame(
            linhas,
            columns=["Iteração", "a", "b", "c", "f(a)", "f(b)", "f(c)", "Erro"],
        ),
        geometria,
    )


def metodo_ponto_fixo(f, g, x0, tol, max_iter):
    linhas = []
    geometria = []
    x = x0

    for i in range(1, max_iter + 1):
        x1 = g(x)
        erro = abs(x1 - x)

        linhas.append([i, x, x1, f(x1), erro])
        geometria.append((x, x1))

        if erro < tol:
            return (
                x1,
                pd.DataFrame(
                    linhas, columns=["Iteração", "xᵢ", "xᵢ₊₁", "f(xᵢ₊₁)", "Erro"]
                ),
                geometria,
            )
        x = x1

    return (
        x,
        pd.DataFrame(
            linhas, columns=["Iteração", "xᵢ", "xᵢ₊₁", "f(xᵢ₊₁)", "Erro"]
        ),
        geometria,
    )


def metodo_newton(f, df, x0, tol, max_iter):
    linhas = []
    geometria = []
    x = x0

    for i in range(1, max_iter + 1):
        fx = f(x)
        dfx = df(x)

        if abs(dfx) < 1e-14:
            raise ValueError("Derivada muito próxima de zero durante Newton-Raphson.")

        x1 = x - fx / dfx
        erro = abs(x1 - x)

        linhas.append([i, x, fx, dfx, x1, erro])
        geometria.append((x, fx, dfx, x1))

        if erro < tol:
            return (
                x1,
                pd.DataFrame(
                    linhas,
                    columns=["Iteração", "xᵢ", "f(xᵢ)", "f'(xᵢ)", "xᵢ₊₁", "Erro"],
