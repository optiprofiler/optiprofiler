(function () {
    var migrationKey = "op_docs_theme_defaulted_to_dark_v1";
    var defaultTheme = "dark";

    function applyDefaultTheme(updateStorage) {
        document.documentElement.dataset.mode = defaultTheme;
        setTheme(defaultTheme, updateStorage);
    }

    function setTheme(theme, updateStorage) {
        document.documentElement.dataset.theme = theme;

        document.querySelectorAll(".dropdown-menu").forEach(function (menu) {
            if (theme === "dark") {
                menu.classList.add("dropdown-menu-dark");
            } else {
                menu.classList.remove("dropdown-menu-dark");
            }
        });

        if (updateStorage) {
            try {
                localStorage.setItem("mode", document.documentElement.dataset.mode || theme);
                localStorage.setItem("theme", theme);
            } catch (error) {
                /* localStorage may be unavailable in private browsing. */
            }
        }

        updateThemeButton(theme);
    }

    function updateThemeButton(theme) {
        document.querySelectorAll(".theme-switch-button").forEach(function (button) {
            var sun = button.querySelector('[data-mode="light"]');
            var moon = button.querySelector('[data-mode="dark"]');
            var auto = button.querySelector('[data-mode="auto"]');

            if (auto) {
                auto.style.display = "none";
            }

            if (sun) {
                sun.style.display = theme === "dark" ? "" : "none";
                sun.title = "Switch to light mode";
            }

            if (moon) {
                moon.style.display = theme === "dark" ? "none" : "";
                moon.title = "Switch to dark mode";
            }

            button.setAttribute(
                "aria-label",
                theme === "dark" ? "Switch to light mode" : "Switch to dark mode"
            );
            button.setAttribute(
                "title",
                theme === "dark" ? "Switch to light mode" : "Switch to dark mode"
            );
            button.dataset.bsTitle =
                theme === "dark" ? "Switch to light mode" : "Switch to dark mode";
        });
    }

    function setManualTheme(theme) {
        document.documentElement.dataset.mode = theme;
        setTheme(theme, true);
    }

    try {
        var storedMode = localStorage.getItem("mode");

        if (storedMode === "light" || storedMode === "dark") {
            document.documentElement.dataset.mode = storedMode;
            setTheme(storedMode, false);
        } else {
            localStorage.setItem(migrationKey, "1");
            applyDefaultTheme(true);
        }
    } catch (error) {
        applyDefaultTheme(false);
    }

    document.addEventListener("DOMContentLoaded", function () {
        updateThemeButton(document.documentElement.dataset.theme || defaultTheme);

        document.querySelectorAll(".theme-switch-button").forEach(function (button) {
            button.addEventListener("click", function (event) {
                event.preventDefault();
                event.stopImmediatePropagation();

                var currentTheme = document.documentElement.dataset.theme || defaultTheme;
                setManualTheme(currentTheme === "dark" ? "light" : "dark");
            }, true);
        });
    });
}());
