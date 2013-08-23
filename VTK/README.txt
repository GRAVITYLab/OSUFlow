OSUFlow-VTK library

This library contains integration classes between VTK and OSUFlow.  The integration includes a new OSUFlowVTK class which allows OSUFlow to load VTK's generic dataset class 'vtkDataSet' perform visualization on top of it, and a new vtkOSUFlow filter which can be integrated in VTK's visualization pipeline.  Each class is saved in a separate file in its name.  Here is the description of each class:

OSUFlowVTK: 
Inherits OSUFlow.  Provides a function to assign vtkDataSet.  

vtkOSUflow:
A VTK filter which inherits vtkStreamer.  Calls OSUFlow by assigning the given seeds and dataset from the input port.  Outputs converted traces in vtkPolyData.

VectorFieldVTK:
An inner class for OSUFlow's flow computation, which inherits and overrides OSUFlow's VectorField class.  Returns vector data by querrying vtkDataSet.



Chun-Ming Chen @ Kitware
May 21, 2013

