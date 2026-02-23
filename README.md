# Atomic orbital modeller

This (proof of concept) site uses WebGL to render hydrogen-like atomic orbitals computed using the three quantum numbers n, l and m. 

Planned improvements:
   - Implementing iso-surfaces
   - Figuring out how to make higher order n not get fainter as the order increases
     - Fixing cut planes
   - Optimising frag.glsl so it runs faster
     - Implementing web workers so the main thread isn't clogged up
   - Snapping camera to axes based on user input
   - Possibly using SQLite to allow for some quality of life settings to be stored semi-persistently.
     - Which of course means exposing those settings to the user in the first place
