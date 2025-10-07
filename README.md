# S4 Slicer
A generic non-planar slicer, that can print almost any part without support.

Please use the [dicussions tab](https://github.com/jyjblrd/S4_Slicer/discussions) to ask questions and help others.

[Try it now](https://colab.research.google.com/github/jyjblrd/S4_Slicer) on Google Colab! (note: colab free tier is only powerful enough to slice very simple models)

[![Watch the video](https://github.com/jyjblrd/S4_Slicer/blob/main/thumnail.jpeg?raw=true)](https://www.youtube.com/watch?v=M51bMMVWbC8)

Check out my [YouTube video](https://youtu.be/M51bMMVWbC8?si=pfud7bHgjYDnO2_z) for more details!

Thank you to JLCCNC for helping create the extruder mount and build plate for my [4 Axis Core R-Theta Printer](https://github.com/jyjblrd/Core-R-Theta-4-Axis-Printer).


## Native C++ path gradient prototype

The original slicer prototype is implemented in Python within `main.ipynb`. For workloads that are
dominated by the `calculate_path_length_to_base_gradient` routine we now ship a standalone, optimised
C++ rewrite in [`cpp/rotation_optimizer.cpp`](cpp/rotation_optimizer.cpp). The implementation keeps the
same inputs as the notebook version (cell centres, neighbourhood lists, face normals, etc.) and performs
multi-source Dijkstra, plane fitting and optional smoothing using cache-friendly data structures.

Build and run the self-contained example with:

```bash
g++ -std=c++17 -O3 cpp/rotation_optimizer.cpp cpp/example.cpp -o path_gradient
./path_gradient
```

The program prints the per-cell gradient, the computed distance to the nearest bottom cell and the
index of that cell, allowing easy verification against the Python workflow. The implementation can be
embedded in other tools or wrapped with Python bindings for tight integration with the existing
notebook pipeline.



Bibtex Citation:
```
@software{Bird_S4_Slicer,
author = {Bird, Joshua},
license = {GPL-3.0},
title = {{S4 Slicer}},
url = {https://github.com/jyjblrd/S4_Slicer}
}
```
