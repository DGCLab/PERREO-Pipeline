#!/usr/bin/env python3
# braille_logo.py
# Uso:
#   python3 braille_logo.py logo.svg 110
#   python3 braille_logo.py logo.png 110
# Opciones:
#   --invert        (invierte tinta/fondo)
#   --edge 1.2      (más contorno, ojos/boca más marcados)
#   --thr 0         (0=adaptativo; o fija un umbral 1..255)
#   --crop 1        (recorte automático al contenido)
#
# Ejemplos:
#   python3 braille_logo.py logo.svg 120 --edge 1.4 --crop 1
#   python3 braille_logo.py logo.svg 90  --edge 1.1 --invert

import argparse
import io
import os
import numpy as np
from PIL import Image

def load_image(path: str, scale: int = 2400) -> Image.Image:
    ext = os.path.splitext(path.lower())[1]
    if ext == ".svg":
        try:
            import cairosvg
        except ImportError as e:
            raise SystemExit("Falta cairosvg. Instala: pip install cairosvg") from e
        # Rasteriza el SVG a PNG grande para mantener detalles
        png_bytes = cairosvg.svg2png(url=path, output_width=scale)
        return Image.open(io.BytesIO(png_bytes)).convert("RGBA")
    else:
        return Image.open(path).convert("RGBA")

def autocrop_rgba(img: Image.Image, pad: int = 6) -> Image.Image:
    arr = np.array(img)
    a = arr[:, :, 3]
    ys, xs = np.where(a > 0)
    if len(xs) == 0:
        return img
    x0, x1 = xs.min(), xs.max()
    y0, y1 = ys.min(), ys.max()
    x0 = max(0, x0 - pad); y0 = max(0, y0 - pad)
    x1 = min(arr.shape[1] - 1, x1 + pad); y1 = min(arr.shape[0] - 1, y1 + pad)
    return img.crop((x0, y0, x1 + 1, y1 + 1))

# Mapeo Braille: cada carácter = 8 puntos (2 columnas x 4 filas)
# Bits (según estándar Braille):
# (0,0)=1, (0,1)=2, (0,2)=4, (0,3)=64
# (1,0)=8, (1,1)=16,(1,2)=32,(1,3)=128
BRAILLE_BASE = 0x2800
BITMAP = np.array([
    [0x01, 0x08],
    [0x02, 0x10],
    [0x04, 0x20],
    [0x40, 0x80],
], dtype=np.uint16)

def to_braille(mask: np.ndarray) -> str:
    """
    mask: array HxW boolean o 0/1 donde 1 = tinta
    Renderiza en Braille con bloques de 4x2 píxeles.
    """
    H, W = mask.shape
    # Asegura múltiplos de 4x2
    H2 = (H // 4) * 4
    W2 = (W // 2) * 2
    mask = mask[:H2, :W2]

    out_lines = []
    for y in range(0, H2, 4):
        line_chars = []
        block = mask[y:y+4, :]  # 4 x W2
        for x in range(0, W2, 2):
            b = block[:, x:x+2].astype(np.uint16)  # 4x2
            code = int((b * BITMAP).sum())
            line_chars.append(chr(BRAILLE_BASE + code))
        out_lines.append("".join(line_chars).rstrip())
    return "\n".join(out_lines)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path", help="logo.svg o logo.png")
    ap.add_argument("width", nargs="?", type=int, default=110, help="ancho en caracteres Braille")
    ap.add_argument("--invert", action="store_true", help="invertir tinta/fondo")
    ap.add_argument("--edge", type=float, default=1.25, help="realce de bordes (1.0-1.6 típico)")
    ap.add_argument("--thr", type=int, default=0, help="0=adaptativo; o umbral fijo 1..255")
    ap.add_argument("--crop", type=int, default=1, help="1=recorte automático a contenido")
    args = ap.parse_args()

    img = load_image(args.path)

    if args.crop:
        img = autocrop_rgba(img)

    # Fondo blanco (evita grises raros por alpha)
    bg = Image.new("RGBA", img.size, (255, 255, 255, 255))
    img = Image.alpha_composite(bg, img).convert("RGB")

    # Resize a ancho en caracteres Braille:
    # Cada char Braille = 2 píxeles en ancho. Por tanto px_width = width*2
    # En alto: cada char = 4 píxeles, así que mantendremos proporción con factor 2/4 = 0.5
    w_px = args.width * 2
    aspect = img.size[1] / img.size[0]
    h_px = max(4, int(w_px * aspect))  # alto en píxeles antes de ajustar a múltiplo de 4
    img = img.resize((w_px, h_px), Image.Resampling.LANCZOS)

    # Procesado: gris + realce contornos + umbral
    import cv2
    gray = cv2.cvtColor(np.array(img), cv2.COLOR_RGB2GRAY)

    # Suaviza sin matar contornos
    gray = cv2.bilateralFilter(gray, d=7, sigmaColor=45, sigmaSpace=45)

    # Unsharp mask suave para ojos/boca/letras
    blur = cv2.GaussianBlur(gray, (0, 0), 1.0)
    sharp = cv2.addWeighted(gray, 1.0 + (args.edge - 1.0), blur, -(args.edge - 1.0), 0)

    # Umbral: adaptativo suele dar texto más legible
    if args.thr and 1 <= args.thr <= 255:
        _, bw = cv2.threshold(sharp, args.thr, 255, cv2.THRESH_BINARY)
    else:
        bw = cv2.adaptiveThreshold(
            sharp, 255,
            cv2.ADAPTIVE_THRESH_GAUSSIAN_C,
            cv2.THRESH_BINARY,
            31, 7
        )

    # Queremos tinta = 1 (oscuro). En THRESH_BINARY, oscuro->0, así que invertimos:
    ink = (bw == 0)

    # Limpieza mínima (quita motas sin engordar)
    ink_u8 = (ink.astype(np.uint8) * 255)
    ink_u8 = cv2.morphologyEx(ink_u8, cv2.MORPH_OPEN, np.ones((2, 2), np.uint8), iterations=1)
    ink = (ink_u8 > 0)

    if args.invert:
        ink = ~ink

    # Asegura múltiplos para Braille
    H, W = ink.shape
    H = (H // 4) * 4
    W = (W // 2) * 2
    ink = ink[:H, :W]

    print(to_braille(ink))

if __name__ == "__main__":
    main()