(function () {
    var migrationKey = "op_docs_theme_defaulted_to_system_v1";

    function getSystemTheme() {
        var prefersDark = window.matchMedia &&
            window.matchMedia("(prefers-color-scheme: dark)").matches;
        return prefersDark ? "dark" : "light";
    }

    function applyAutoTheme(updateStorage) {
        var theme = getSystemTheme();

        document.documentElement.dataset.mode = "auto";
        setTheme(theme, updateStorage);
        document.documentElement.dataset.mode = "auto";
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
        var migrated = localStorage.getItem(migrationKey);

        if (!migrated) {
            localStorage.setItem(migrationKey, "1");
            applyAutoTheme(true);
        } else if (!storedMode || storedMode === "auto") {
            applyAutoTheme(!storedMode);
        }
    } catch (error) {
        applyAutoTheme(false);
    }

    if (window.matchMedia) {
        var query = window.matchMedia("(prefers-color-scheme: dark)");
        var updateAutoMode = function () {
            try {
                var mode = localStorage.getItem("mode");
                if (mode && mode !== "auto") {
                    return;
                }
            } catch (error) {
                /* Fall through and update the DOM. */
            }
            applyAutoTheme(true);
        };

        if (query.addEventListener) {
            query.addEventListener("change", updateAutoMode);
        } else if (query.addListener) {
            query.addListener(updateAutoMode);
        }
    }

    document.addEventListener("DOMContentLoaded", function () {
        updateThemeButton(document.documentElement.dataset.theme || getSystemTheme());

        document.querySelectorAll(".theme-switch-button").forEach(function (button) {
            button.addEventListener("click", function (event) {
                event.preventDefault();
                event.stopImmediatePropagation();

                var currentTheme = document.documentElement.dataset.theme || getSystemTheme();
                setManualTheme(currentTheme === "dark" ? "light" : "dark");
            }, true);
        });
    });
}());
