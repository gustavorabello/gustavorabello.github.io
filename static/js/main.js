document.addEventListener("DOMContentLoaded", () => {
  const navButton = document.querySelector("[data-nav-toggle]");
  const navMenu = document.querySelector("[data-nav-menu]");

  if (navButton && navMenu) {
    navButton.addEventListener("click", () => {
      navMenu.classList.toggle("is-open");
    });
  }

  // In-page anchor links (e.g. the table of contents and citations) use bare
  // "#id" fragments. Because every page sets a <base href> so relative asset
  // paths resolve at any depth, the browser resolves those fragments against
  // the base URL (the site root) and clicking them jumps to the homepage.
  // Rewrite them to point at the current page so they scroll in-page instead.
  if (document.querySelector("base")) {
    const here = window.location.pathname + window.location.search;
    document.querySelectorAll('a[href^="#"]').forEach((link) => {
      const fragment = link.getAttribute("href");
      if (fragment.length > 1) {
        link.setAttribute("href", here + fragment);
      }
    });
  }
});
