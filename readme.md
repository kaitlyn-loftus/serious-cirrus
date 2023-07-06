## Idealized 3D GCM and microphysical modeling of water ice clouds on early Mars: Implications for climate 

Code for microphysics calculations and their visualizations associated with "Idealized 3D GCM and microphysical modeling of water ice clouds on early Mars: Implications for climate" by Ding, Loftus, & Wordsworth (submitted). See [FMSPCM](https://github.com/fdingdfdfdf/FMSPCM) for GCM code component.

Steps to reproduce microhpysics calculations & visualizations:

1. download repo

2. in repo directory, start Julia (v1.9 used originally)

3. in Julia REPL:
```
] activate .
instantiate
```
   (backspace to exit Pkg REPL)
   ```
   include(“calc.jl”)
   include(“plot.jl”)
   ```

See ``src/MΦ.jl`` for the majority of calculation details.
  

Please contact [Kaitlyn Loftus](mailto:kaitlyn.loftus@columbia.edu) with questions / issues / concerns.