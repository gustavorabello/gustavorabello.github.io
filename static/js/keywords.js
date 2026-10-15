(() => {
  const root = document.getElementById("research-topics");
  if (!root) return;
  const dataElement = document.getElementById("keyword-data");
  if (!dataElement) return;
  const data = JSON.parse(dataElement.textContent);
  if (!data.terms.length) return;
  const cloud = root.querySelector(".keyword-cloud");
  const words = [...cloud.querySelectorAll(".keyword-word")];
  const results = root.querySelector(".keyword-results");
  const heading = document.getElementById("keyword-result-title");
  const filters = root.querySelector(".keyword-categories");
  const timeline = root.querySelector(".keyword-timeline");
  const resultLayout = root.querySelector(".keyword-results-layout");
  const videoSection = root.querySelector(".keyword-videos");
  const videoGrid = root.querySelector(".keyword-video-grid");
  const reducedMotion = matchMedia("(prefers-reduced-motion: reduce)");
  let selected = null;
  let category = "all";
  let categoryRevision = 0;
  let viewTransition = null;

  function node(tag, className, text) {
    const element = document.createElement(tag);
    element.className = className;
    if (text != null) element.textContent = text;
    return element;
  }

  function saveSelection() {
    const params = new URLSearchParams({topic: selected.id});
    if (category !== "all") params.set("type", category);
    history.replaceState(null, "", `${location.pathname}${location.search}#${params}`);
  }

  function renderTimeline() {
    const videoView = category === "videos";
    resultLayout.classList.toggle("is-video-view", videoView);
    timeline.hidden = videoView;
    filters.querySelectorAll("button").forEach(button => {
      button.setAttribute("aria-pressed", String(button.dataset.category === category));
    });
    if (videoView) {
      const videos = selected.video_items.map(index => data.videos[index]);
      const years = new Set(videos.map(video => video.year));
      root.querySelector(".keyword-total").textContent = `${videos.length} related ${videos.length === 1 ? "video" : "videos"} · ${years.size} ${years.size === 1 ? "year" : "years"}`;
      return;
    }
    const items = selected.items.map(index => data.items[index])
      .filter(item => category === "all" || item.category === category)
      .sort((a, b) => Number(b.year) - Number(a.year) || a.category.localeCompare(b.category) || a.title.localeCompare(b.title));
    const fragment = document.createDocumentFragment();
    const years = [...new Set(items.map(item => item.year))];
    for (const year of years) {
      const entries = items.filter(item => item.year === year);
      const group = node("section", "keyword-year");
      group.dataset.year = year;
      group.setAttribute("aria-label", `${year}: ${entries.length} related items`);
      const label = node("h4", "keyword-year-heading", year);
      label.append(node("span", "keyword-year-count", `${entries.length} ${entries.length === 1 ? "item" : "items"}`));
      const cards = node("div", "keyword-year-cards");
      for (const item of entries) {
        const card = node("a", `keyword-card ${item.kind_class}`);
        card.href = item.url;
        card.setAttribute("aria-label", `${item.type_label}: ${item.title}`);
        if (category === "all") card.append(node("p", "keyword-card-type", item.type_label));
        card.append(node("h5", "keyword-card-title", item.title));
        if (item.venue) card.append(node("p", "keyword-card-venue", item.venue));
        cards.append(card);
      }
      group.append(label, cards);
      fragment.append(group);
    }
    timeline.replaceChildren(fragment);
    root.querySelector(".keyword-total").textContent = category === "all"
      ? `${selected.count} related items · ${years.length} years`
      : `${items.length} of ${selected.count} related items · ${years.length} years`;
  }

  function changeCategory(nextCategory) {
    if (nextCategory === category) return;
    const revision = ++categoryRevision;
    const topic = selected;
    const changingVideoView = (nextCategory === "videos") !== (category === "videos");
    viewTransition?.skipTransition();
    const update = () => {
      if (revision !== categoryRevision || selected !== topic) return;
      category = nextCategory;
      renderTimeline();
      saveSelection();
      if (nextCategory === "videos") results.scrollIntoView({behavior: "instant", block: "start"});
    };
    if (changingVideoView && !reducedMotion.matches && document.startViewTransition) {
      viewTransition = document.startViewTransition(update);
    } else {
      update();
      if (changingVideoView && !reducedMotion.matches) {
        videoGrid.animate([{opacity: .35, transform: "translateX(24px)"}, {opacity: 1, transform: "translateX(0)"}],
          {duration: 300, easing: "ease-out"});
      }
    }
  }

  function renderVideos() {
    const videos = (selected.video_items || []).map(index => data.videos[index]);
    videoSection.hidden = !videos.length;
    resultLayout.classList.toggle("has-videos", videos.length > 0);
    root.querySelector(".keyword-video-total").textContent = `${videos.length} ${videos.length === 1 ? "video" : "videos"}`;
    const fragment = document.createDocumentFragment();
    for (const [index, video] of videos.entries()) {
      const card = node("a", "keyword-video-card");
      card.style.viewTransitionName = `keyword-video-${index}`;
      card.href = video.url;
      card.title = video.title;
      card.setAttribute("aria-label", `Watch video: ${video.title}`);
      const image = node("img", "keyword-video-thumbnail");
      image.src = video.thumbnail;
      image.alt = "";
      image.width = 160;
      image.height = 90;
      image.loading = "lazy";
      const copy = node("div", "keyword-video-copy");
      copy.append(node("p", "keyword-video-meta", `${video.year}${video.duration ? ` · ${video.duration}` : ""}`));
      copy.append(node("h5", "keyword-video-title", video.title));
      card.append(image, copy);
      fragment.append(card);
    }
    videoGrid.replaceChildren(fragment);
    videoGrid.scrollTop = 0;
  }

  function selectTopic(term, requestedCategory = "all", scroll = false) {
    categoryRevision++;
    viewTransition?.skipTransition();
    selected = term;
    const available = data.categories.filter(group => term.items.some(index => data.items[index].category === group.id));
    available.push({id: "videos", label: "Videos", kind_class: "is-video"});
    const videoCount = (term.video_items || []).length;
    category = available.some(group => group.id === requestedCategory) && (requestedCategory !== "videos" || videoCount)
      ? requestedCategory : "all";
    words.forEach(word => word.setAttribute("aria-pressed", String(word.dataset.keyword === term.id)));
    heading.textContent = term.label;
    filters.replaceChildren();
    for (const group of [{id: "all", label: "All types", kind_class: ""}, ...available]) {
      const count = group.id === "all" ? term.count : group.id === "videos" ? videoCount
        : term.items.filter(index => data.items[index].category === group.id).length;
      const button = node("button", `keyword-category ${group.kind_class}`, group.label);
      button.type = "button";
      button.dataset.category = group.id;
      button.disabled = group.id === "videos" && !videoCount;
      button.append(node("span", "keyword-category-count", String(count)));
      button.addEventListener("click", () => changeCategory(group.id));
      filters.append(button);
    }
    renderTimeline();
    renderVideos();
    results.hidden = false;
    if (scroll) {
      heading.focus({preventScroll: true});
      results.scrollIntoView({behavior: reducedMotion.matches ? "instant" : "smooth", block: "start"});
    }
  }

  words.forEach(word => {
    const term = data.terms.find(item => item.id === word.dataset.keyword);
    word.title = `${term.count} projects and publications · ${term.video_count || 0} videos`;
    word.addEventListener("click", () => {
      selectTopic(term, "all", true);
      saveSelection();
    });
  });

  function restoreSelection() {
    const params = new URLSearchParams(location.hash.slice(1));
    const term = data.terms.find(item => item.id === params.get("topic"));
    if (term) selectTopic(term, params.get("type") || "all");
    else {
      results.hidden = true;
      words.forEach(word => word.setAttribute("aria-pressed", "false"));
    }
  }
  restoreSelection();
  window.addEventListener("hashchange", restoreSelection);

  let previousWidth = 0;
  const textMeasure = document.createElement("canvas").getContext("2d");
  function layoutCloud() {
    const width = cloud.clientWidth;
    if (!width) return;
    cloud.classList.add("is-positioned");
    const coarse = matchMedia("(pointer: coarse)").matches;
    const maximum = Math.sqrt(data.terms[0].count);
    const minimum = Math.sqrt(data.terms[data.terms.length - 1].count);
    const integrated = root.classList.contains("keyword-explorer-integrated");
    let height = integrated ? (coarse && width < 540 ? 360 : width < 400 ? 280 : 220) : width < 540 ? 500 : 420;
    // Retry with denser text or a taller field rather than hiding smaller topics.
    for (let attempt = 0; attempt < 12; attempt++) {
      const placed = [];
      const scale = Math.max(.7, 1 - attempt * .05);
      cloud.style.height = `${height}px`;
      let complete = true;
      for (let index = 0; index < words.length; index++) {
        const word = words[index];
        word.style.left = "0px";
        word.style.top = "0px";
        const weight = maximum === minimum ? 0 : (Math.sqrt(data.terms[index].count) - minimum) / (maximum - minimum);
        const leadingSize = integrated ? Math.min(44, width / 10) : Math.min(62, width / 10);
        const surroundingSize = integrated ? 11 + weight * (width < 400 ? 10 : 16) : 13 + weight * (width < 540 ? 10 : 32);
        let size = index === 0 ? leadingSize : Math.max(11, surroundingSize * scale);
        const style = getComputedStyle(word);
        const measure = () => {
          textMeasure.font = `${style.fontWeight} ${size}px ${style.fontFamily}`;
          return Math.ceil(textMeasure.measureText(word.textContent).width - .04 * size * Math.max(0, word.textContent.length - 1) + 8);
        };
        let wordWidth = measure();
        if (wordWidth > width - 32) {
          size *= (width - 32) / wordWidth;
          wordWidth = measure();
        }
        word.style.fontSize = `${size}px`;
        const box = {width: wordWidth, height: Math.max(Math.ceil(size * 1.05 + 6), coarse ? 44 : 0)};
        word.style.width = `${box.width}px`;
        word.style.height = `${box.height}px`;
        const w = box.width + 8;
        const h = Math.max(box.height + 8, coarse ? 48 : 0);
        let found = false;
        for (let step = 0; step < 8000; step++) {
          const radius = index === 0 ? 0 : 4.7 * Math.sqrt(step);
          const angle = step * .39 + index * 1.9;
          const x = width / 2 + radius * Math.cos(angle) * width / height - w / 2;
          const y = height / 2 + radius * Math.sin(angle) - h / 2;
          if (x < 7 || y < 7 || x + w > width - 7 || y + h > height - 7) continue;
          if (placed.some(rect => x < rect.x + rect.w && x + w > rect.x && y < rect.y + rect.h && y + h > rect.y)) continue;
          word.style.left = `${x + 4}px`;
          word.style.top = `${y + (h - box.height) / 2}px`;
          placed.push({x, y, w, h});
          found = true;
          break;
        }
        if (!found) { complete = false; break; }
      }
      if (complete) return;
      if (attempt >= 3) height += 70;
    }
    // Keep every keyword accessible if the font/viewport defeats packing.
    cloud.classList.remove("is-positioned");
    cloud.style.height = "auto";
  }
  const observer = new ResizeObserver(entries => {
    const width = Math.round(entries[0].contentRect.width);
    if (width !== previousWidth) { previousWidth = width; layoutCloud(); }
  });
  observer.observe(cloud);
  layoutCloud();
  if (document.fonts) document.fonts.ready.then(layoutCloud);
})();
