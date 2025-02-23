from shiny import reactive, render
from shiny.express import input, ui

import asyncio
# --------------------------------------------------------------------
ui.page_opts(title = "📡 some title here",font="bold")
# --------------------------------------------------------------------
with ui.sidebar(open="desktop",fillable=True):
  ui.input_action_button("load_data", "load data", font="bold")
  ui.input_action_button("compute_now", "⟶ compute now", font="bold")
  ui.input_action_button("reset", "reset", font="bold")
# --------------------------------------------------------------------
with ui.layout_columns(col_widths=[6, 6, 12]):
  with ui.card(full_screen=True):
    ui.card_header("some plot here",font="bold")
    @render.plot
    def plot():
      import matplotlib.pyplot as plt
      return plt.scatter([1, 2, 3], [4, 5, 6])
# --------------------------------------------------------------------
  with ui.card(full_screen=True):
    ui.card_header("some other plot here too",font="bold")
    @render.plot
    @reactive.event(input.compute_now)
    def ploto():
      import matplotlib.pyplot as plt
      return plt.scatter([4, 5, 6], [1, 2, 3])
# --------------------------------------------------------------------
# 
# 
# 
# --------------------------------------------------------------------
@render.ui
@reactive.event(input.compute_now)
async def compute():
  with ui.Progress(min=1, max=15) as p:
    p.set(message="computing...", detail="")

    for i in range(1, 15):
      p.set(i, message="computing...")
      await asyncio.sleep(0.1)

  return " :: computing done ::"
# --------------------------------------------------------------------
# --------------------------------------------------------------------