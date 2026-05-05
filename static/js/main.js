document.addEventListener("DOMContentLoaded", () => {
  const navButton = document.querySelector("[data-nav-toggle]");
  const navMenu = document.querySelector("[data-nav-menu]");

  if (navButton && navMenu) {
    navButton.addEventListener("click", () => {
      navMenu.classList.toggle("is-open");
    });
  }
});
