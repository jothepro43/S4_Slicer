# S4 Slicer
A generic non-planar slicer, that can print almost any part without support.

Please use the [dicussions tab](https://github.com/jyjblrd/S4_Slicer/discussions) to ask questions and help others.

[Try it now](https://colab.research.google.com/github/jyjblrd/S4_Slicer) on Google Colab! (note: colab free tier is only powerful enough to slice very simple models)

[![Watch the video](https://github.com/jyjblrd/S4_Slicer/blob/main/thumnail.jpeg?raw=true)](https://www.youtube.com/watch?v=M51bMMVWbC8)

Check out my [YouTube video](https://youtu.be/M51bMMVWbC8?si=pfud7bHgjYDnO2_z) for more details!

Thank you to JLCCNC for helping create the extruder mount and build plate for my [4 Axis Core R-Theta Printer](https://github.com/jyjblrd/Core-R-Theta-4-Axis-Printer).


## Native C++ mesh processing toolkit

The original slicer prototype is implemented in Python within `main.ipynb`. For workloads that are
dominated by the `calculate_path_length_to_base_gradient` routine we now ship a standalone, optimised
C++ rewrite in [`cpp/rotation_optimizer.cpp`](cpp/rotation_optimizer.cpp). The implementation keeps the
same inputs as the notebook version (cell centres, neighbourhood lists, face normals, etc.) and performs
multi-source Dijkstra, plane fitting and optional multi-threaded smoothing using cache-friendly data
structures.

Two small utilities demonstrate how to embed the native routines:

* [`cpp/example.cpp`](cpp/example.cpp) keeps the original, in-memory smoke test which can be compiled with:

  ```bash
  g++ -std=c++17 -O3 cpp/rotation_optimizer.cpp cpp/example.cpp -o path_gradient
  ./path_gradient
  ```

* [`cpp/mesh_processor.cpp`](cpp/mesh_processor.cpp) loads an arbitrary STL mesh, constructs the required
  connectivity on the fly and exports both a cleaned STL and a coloured PLY heat-map visualisation of the
  gradient. The build requires the STL utilities in [`cpp/stl_mesh.cpp`](cpp/stl_mesh.cpp):

  ```bash
  g++ -std=c++17 -O3 -pthread cpp/mesh_processor.cpp cpp/stl_mesh.cpp cpp/rotation_optimizer.cpp -o mesh_processor
  ./mesh_processor input.stl processed.stl gradient.ply
  ```

  The optional final argument specifies the maximum overhang angle in degrees (defaults to 30°). The
  generated PLY can be inspected in MeshLab, Blender or any viewer that supports vertex colours.

Both binaries are ignored by git and may be safely rebuilt locally whenever the underlying algorithms
change.



Bibtex Citation:
```
@software{Bird_S4_Slicer,
author = {Bird, Joshua},
license = {GPL-3.0},
title = {{S4 Slicer}},
url = {https://github.com/jyjblrd/S4_Slicer}
}
```
