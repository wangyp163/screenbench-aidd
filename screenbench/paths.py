"""工具链解析：vina / obabel 用本地 deb 解包版（tools/local），无需 root。"""
import os
import shutil
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
TOOLS = REPO / "tools" / "local" / "usr"


def tool_bin(name: str) -> str:
    """返回工具绝对路径；优先 tools/local，其次 PATH，最后报可执行错误。"""
    local = TOOLS / "bin" / name
    if local.is_file() and os.access(local, os.X_OK):
        return str(local)
    found = shutil.which(name)
    if found:
        return found
    raise RuntimeError(f"找不到可执行文件 {name}；请先运行 tools 下的 deb 解包（见 PITFALLS.md）")


def tool_env() -> dict:
    """vina/obabel 运行所需环境（LD_LIBRARY_PATH + BABEL_LIBDIR 指向本地解包库）。"""
    env = os.environ.copy()
    lib = TOOLS / "lib" / "x86_64-linux-gnu"
    if lib.is_dir():
        env["LD_LIBRARY_PATH"] = f"{lib}:{env.get('LD_LIBRARY_PATH', '')}"
    babel = lib / "openbabel" / "3.1.1"
    if babel.is_dir():
        env["BABEL_LIBDIR"] = str(babel)
    return env
