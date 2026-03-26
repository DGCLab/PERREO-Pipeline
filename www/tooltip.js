function initPerreoTooltips() {
  const triggers = [].slice.call(
    document.querySelectorAll('[data-bs-toggle="tooltip"]')
  );

  triggers.forEach((el) => {
    if (!bootstrap.Tooltip.getInstance(el)) {
      new bootstrap.Tooltip(el);
    }
  });
}

document.addEventListener('DOMContentLoaded', initPerreoTooltips);
document.addEventListener('shown.bs.tab', initPerreoTooltips);
document.addEventListener('shiny:inputchanged', initPerreoTooltips);
