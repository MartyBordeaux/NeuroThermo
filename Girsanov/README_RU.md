# Girsanov for NeuroThermo: complete available-support analysis

Произведён модельный расчёт для всех 264 frozen retained WT×SCA3 сценариев; итоговый научный документ `SCIENTIFIC_REPORT_RU.md`.

Каталоги:
- `originals/` — неизменённые сохранённые входы: 34 retained параметризации, архивные pair-specific маркеры диагонального пути и прежний 31-точечный coupled расчет.
- `all264_2d_D/` — 11×11 тензор Фишера, 11 известных предписанных маршрутов × 264 сценария, 4 независимые модели в точке.
- `endpoint_all264/` — точный ожидаемый finite-horizon Girsanov KL между endpoint-законами для 5 условий протокола/шума × 2 reference-направления × 264 сценария; 16 независимых траекторий на law.
- `endpoint_precision64/` — N=64 sensitivity для 17 отобранных сценариев с высокой MC вариабельностью.
- `trace_qc/` — реальные x(t),y(t),z(t), spike count и картинки для трех репрезентативных best/best пар.
- `final_analysis/` — pair-weighted таблицы и графики, сравнение диагонали с frozen firing и numerical QC.
- `full2d.py`, `endpoint_kl.py`, `summarize.py`, `quality_tests.py`, `trace_qc.py`, `refine_highcv.py` — исходники.

Команды, идентичные результатам, см. в конце SCIENTIFIC_REPORT_RU.md. Результаты нельзя трактовать как проверку оригинальных 525 frozen маршрутов: исходный маршрутный бинарник SHA256 132d85c9fe10c5b9e50e2cbf0e9ebcb816c084deadf88e8841ec31f26e470d5e отсутствует.
