document.addEventListener("DOMContentLoaded", () => {
  const navButton = document.querySelector("[data-nav-toggle]");
  const navMenu = document.querySelector("[data-nav-menu]");

  if (navButton && navMenu) {
    navButton.addEventListener("click", () => {
      navMenu.classList.toggle("is-open");
      const isOpen = navMenu.classList.contains("is-open");
      navButton.setAttribute("aria-expanded", isOpen);
      navButton.setAttribute("aria-label", isOpen ? "Close navigation" : "Open navigation");
    });
  }

  const revealItems = document.querySelectorAll("[data-reveal]");
  if ("IntersectionObserver" in window && !window.matchMedia("(prefers-reduced-motion: reduce)").matches) {
    document.documentElement.classList.add("has-reveal");
    const revealObserver = new IntersectionObserver((entries) => {
      entries.forEach((entry) => {
        if (entry.isIntersecting) {
          entry.target.classList.add("is-visible");
          revealObserver.unobserve(entry.target);
        }
      });
    }, { rootMargin: "0px 0px -7%", threshold: 0.01 });
    revealItems.forEach((item) => revealObserver.observe(item));
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
