(function() {
    const footer = document.createElement('footer');
    footer.className = 'footer';
    footer.innerHTML = '© 2026 Alexander Sandercock. Licensed under <a href="https://www.gnu.org/licenses/gpl-3.0.html" target="_blank" rel="noopener">GPL-3.0</a>.';

    const content = document.querySelector('.content');
    if (content) {
        content.appendChild(footer);
    }
})();
