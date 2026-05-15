(function () {
    function applySystemTheme() {
        var prefersDark = window.matchMedia &&
            window.matchMedia("(prefers-color-scheme: dark)").matches;
        var theme = prefersDark ? "dark" : "light";

        document.documentElement.dataset.mode = "auto";
        document.documentElement.dataset.theme = theme;

        try {
            localStorage.setItem("mode", "auto");
            localStorage.setItem("theme", theme);
        } catch (error) {
            /* localStorage may be unavailable in private browsing. */
        }
    }

    applySystemTheme();

    if (window.matchMedia) {
        var query = window.matchMedia("(prefers-color-scheme: dark)");
        if (query.addEventListener) {
            query.addEventListener("change", applySystemTheme);
        } else if (query.addListener) {
            query.addListener(applySystemTheme);
        }
    }

    document.addEventListener("click", function (event) {
        if (event.target.closest(".theme-switch-button")) {
            window.setTimeout(applySystemTheme, 0);
        }
    }, true);
}());
