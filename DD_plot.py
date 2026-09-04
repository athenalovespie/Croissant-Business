"""Plot test_DD_psi finite-difference errors with a log-log best-fit line."""

import csv
import math
from pathlib import Path


project_folder = Path(__file__).resolve().parent
data_file = project_folder / "DD_plot_data.csv"
image_file = project_folder / "DD_plot.svg"

x_axis = []
y_axis = []

with data_file.open(newline="", encoding="utf-8") as file:
    for row in csv.DictReader(file):
        x = float(row["x_axis"])
        y = float(row["y_axis"])
        if x > 0.0 and y > 0.0:
            x_axis.append(x)
            y_axis.append(y)

if not x_axis:
    raise ValueError("DD_plot_data.csv does not contain positive plotting points.")

log_x = [math.log10(value) for value in x_axis]
log_y = [math.log10(value) for value in y_axis]
mean_x = sum(log_x) / len(log_x)
mean_y = sum(log_y) / len(log_y)
x_variance = sum((x - mean_x) ** 2 for x in log_x)
if x_variance == 0.0:
    raise ValueError("A best-fit line requires at least two different x-values.")

slope = sum((x - mean_x) * (y - mean_y) for x, y in zip(log_x, log_y)) / x_variance
intercept = mean_y - slope * mean_x

width, height = 900, 600
left, right, top, bottom = 100, 40, 60, 90
plot_width = width - left - right
plot_height = height - top - bottom
x_min, x_max = min(log_x), max(log_x)
fit_y_at_min = slope * x_min + intercept
fit_y_at_max = slope * x_max + intercept
y_min = min(log_y + [fit_y_at_min, fit_y_at_max])
y_max = max(log_y + [fit_y_at_min, fit_y_at_max])


def screen_x(value):
    return left + (value - x_min) / (x_max - x_min) * plot_width


def screen_y(value):
    return top + (y_max - value) / (y_max - y_min) * plot_height


svg = [
    f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
    '<rect width="100%" height="100%" fill="white"/>',
    '<text x="450" y="32" text-anchor="middle" font-family="sans-serif" font-size="22">DD_psi Finite-Difference Error</text>',
]

for tick in range(math.floor(x_min), math.ceil(x_max) + 1):
    if x_min <= tick <= x_max:
        x = screen_x(tick)
        svg.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{height-bottom}" stroke="#dddddd"/>')
        svg.append(f'<text x="{x:.2f}" y="{height-bottom+25}" text-anchor="middle" font-family="sans-serif">10^{tick}</text>')

for tick in range(math.floor(y_min), math.ceil(y_max) + 1):
    if y_min <= tick <= y_max:
        y = screen_y(tick)
        svg.append(f'<line x1="{left}" y1="{y:.2f}" x2="{width-right}" y2="{y:.2f}" stroke="#dddddd"/>')
        svg.append(f'<text x="{left-12}" y="{y+5:.2f}" text-anchor="end" font-family="sans-serif">10^{tick}</text>')

svg.extend([
    f'<line x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}" stroke="black" stroke-width="2"/>',
    f'<line x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}" stroke="black" stroke-width="2"/>',
    f'<text x="{left + plot_width/2}" y="{height-30}" text-anchor="middle" font-family="sans-serif">Epsilon (x_axis)</text>',
    f'<text x="25" y="{top + plot_height/2}" text-anchor="middle" font-family="sans-serif" transform="rotate(-90 25 {top + plot_height/2})">Frobenius error (y_axis)</text>',
])

for x, y in zip(log_x, log_y):
    svg.append(f'<circle cx="{screen_x(x):.2f}" cy="{screen_y(y):.2f}" r="4" fill="#1665d8" opacity="0.8"/>')

svg.append(
    f'<line x1="{screen_x(x_min):.2f}" y1="{screen_y(fit_y_at_min):.2f}" '
    f'x2="{screen_x(x_max):.2f}" y2="{screen_y(fit_y_at_max):.2f}" '
    'stroke="#e34a33" stroke-width="3"/>'
)
svg.extend([
    '<rect x="115" y="72" width="260" height="58" fill="white" stroke="#bbbbbb"/>',
    '<circle cx="132" cy="91" r="4" fill="#1665d8"/>',
    '<text x="145" y="96" font-family="sans-serif">DD_psi error</text>',
    '<line x1="122" y1="115" x2="142" y2="115" stroke="#e34a33" stroke-width="3"/>',
    f'<text x="150" y="120" font-family="sans-serif">Best fit (slope = {slope:.4f})</text>',
    '</svg>',
])

image_file.write_text("\n".join(svg), encoding="utf-8")
print(f"Plotted {len(x_axis)} DD_psi points.")
print(f"DD_psi best-fit slope in log-log space: {slope:.6f}")
print(f"Saved plot to: {image_file}")
