/**
@page step_54 The step-54 tutorial program
This tutorial depends on step-53.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#CADsurfaces"> CAD surfaces </a>
        <li><a href="#TheCADboundaryprojectorclasses"> The CAD boundary projector classes </a>
        <li><a href="#Thetestcase"> The testcase </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#TheTriangulationOnCADclass">The TriangulationOnCAD class</a>
      <ul>
        <li><a href="#TriangulationOnCADTriangulationOnCAD">TriangulationOnCAD::TriangulationOnCAD</a>
        <li><a href="#TriangulationOnCADread_domain">TriangulationOnCAD::read_domain</a>
        <li><a href="#TriangulationOnCADrefine_mesh">TriangulationOnCAD::refine_mesh</a>
        <li><a href="#TriangulationOnCADoutput_results">TriangulationOnCAD::output_results</a>
        <li><a href="#TriangulationOnCADrun">TriangulationOnCAD::run</a>
      </ul>
        <li><a href="#Themainfunction">The main() function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-54/doc/intro.dox

 <br> 

<i>This program was contributed by Andrea Mola and Luca Heltai.</i>

 @note  这个程序阐述了工业几何的概念，使用与OpenCASCADE库（http://www.opencascade.org）接口的工具，允许指定任意的IGES文件来描述你的几何形状的边界。

 @dealiiTutorialDOI{10.5281/zenodo.546220,https://zenodo.org/badge/DOI/10.5281/zenodo.546220.svg} 

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>



在之前的一些教程中（第1步、第3步、第5步、第6步和第49步等），我们已经学会了如何使用deal.II中提供的网格细化方法。这些教程展示了如何利用这些工具为一次模拟产生一个精细的网格，如步骤3；或者从一个粗大的网格开始，在自适应细化的网格上进行一系列模拟，如步骤6的情况。无论采取哪种方法，网格细化都需要对计算域边界进行适当的几何描述，以便在每次细化时将新的网格节点放到边界面上。例如，第5步显示了如何创建一个圆形网格，将一个圆形流形对象自动附加到计算域上，从而使位于边界上的面被细化到圆形上。第53步显示了如何用一个由实验获得的数据定义的流形来做这件事。但是，至少就基本边界形状而言，deal.II实际上只提供了圆、球、盒和其他基本组合。在本教程中，我们将展示如何使用一组开发的类来导入任意的CAD几何图形，将它们分配到计算域的所需边界，并在这种复杂的形状上细化计算网格。




<a name="CADsurfaces"></a><h3> CAD surfaces </h3>


在最常见的工业实践中，任意形状的物体的几何模型是通过计算机辅助设计（CAD）工具实现的。在过去的几十年里，CAD建模器的使用已经普及，因为它们可以为每个设计对象生成一个完整的虚拟模型，通过计算机可以在实物制作之前对其最精细的细节进行可视化、检查和分析。  从数学的角度来看，CAD建模人员的引擎是由分析几何学来代表的，特别是由参数化的曲线和曲面，如B-splines和NURBS，它们足够丰富，可以代表大多数的实际利益的表面。  一旦一个虚拟模型准备好了，所需物体的所有几何特征都被存储在文件中，这些文件实质上包含了构成该物体的参数化曲面和曲线的系数。根据用于定义几何模型的具体CAD工具，当然有几种不同的文件格式，可以组织CAD模型的信息。为了提供一个跨CAD工具交换数据的共同基础，美国国家标准局在1980年发布了初始图形交换代表（IGES）中性文件格式，在本例中使用。

<a name="TheCADboundaryprojectorclasses"></a><h3> The CAD boundary projector classes </h3>


为了导入和查询CAD模型，deal.II库为CAD建模的OpenCASCADE开源库实现了一系列的包装函数。这些函数允许将IGES文件导入OpenCASCADE本地对象，并将其包裹在一系列Manifold类中。

一旦从IGES文件导入，模型就被存储在一个 <code>TopoDS_Shape</code> 中，这是OpenCASCADE框架中定义的通用拓扑实体。从 <code>TopoDS_Shape</code> 中，就可以访问构成它的所有子形状（如顶点、边和面），以及它们的几何描述。在deal.II框架中，组成一个形状的拓扑实体被用来创建一个相应的Manifold表示。在步骤6中，我们看到了如何使用 GridGenerator::hyper_sphere() 来创建一个超球体，它自动将一个SphericalManifold附加到所有边界面。这保证了边界面在网格细化过程中保持在球体或圆上。CAD建模界面的功能被设计为保留相同的结构，允许用户使用导入的CAD形状建立一个投影仪对象，保持我们在其他教程程序中使用的相同程序，即把这种投影仪对象分配给粗略网格的单元、面或边。在每个细化周期，新的网格节点将通过将现有对象的中点投影到指定的几何体上而自动生成。

与球形或圆形边界不同，具有复杂几何形状的边界带来的问题是，在规定形状上细化后创建的新节点最好放在哪里。例如，PolarManifold将周围的点转换为极坐标，计算该坐标系中的平均值（对每个坐标单独计算），最后将点转换回直角坐标。

不过，在一个任意的复杂形状的情况下，一个合适的新节点的位置选择不可能那么容易确定。deal.II中的OpenCASCADE封装器提供了几个采用不同投影策略的投影仪类。第一个投影仪，在 OpenCASCADE::ArclengthProjectionLineManifold 类中实现，只用于边缘细化。它的建立是给它分配一个维度为1的拓扑形状，或者是一个 <code>TopoDS_Edge</code> or a <code>TopoDS_Wire</code> （这是一个复合形状，由几个连接的 <code>TopoDS_Edge</code> 组成），并细化网格边缘，找到新的顶点作为点，将CAD曲线部分的曲线长度分成两个偶数部分，位于原始边缘的顶点之间。

 <img src="https://www.dealii.org/images/steps/developer/step-54.CurveSplit.png" alt="" width="500"> 


在 OpenCASCADE::NormalProjectionBoundary 类中实现了一个不同的投影策略。在构造时分配的 <code>TopoDS_Shape</code> 可以是任意的（图形、面、边的集合或单个面或边都可以）。新的单元格节点首先通过对周围的点进行平均计算，方法与FlatManifold相同。在第二步中，所有的新节点将沿着形状的法线方向被投射到 <code>TopoDS_Shape</code> 。如果没有法线投影，则选择最接近形状的点--通常位于形状的边界上--。  如果形状是由几个子形状组成的，则投影到每个子形状上，并选择最近的投影点。

 <img src="https://www.dealii.org/images/steps/developer/step-54.NormalProjectionEdge.png" alt="" width="500">  <img src="https://www.dealii.org/images/steps/developer/step-54.NormalProjection.png" alt="" width="500">  。

正如我们即将体验到的，对于某些形状，将投影方向设置为CAD表面的法线，将不会导致合适质量的表面网格元素。这是因为CAD表面的法线方向原则上与网格需要新节点所在的方向无关。在这种情况下， OpenCASCADE::DirectionalProjectionBoundary 类可以提供帮助。这个类的构造是指定一个 <code>TopoDS_Shape</code> （至少包含一个面）和一个方向，所有的投影将沿着这个方向进行。新的点将被计算出来，首先对周围的点进行平均化（就像FlatManifold的情况一样），然后沿着构造时使用的方向，在拓扑形状和通过所得到的点的线之间取得最近的交点。  这样一来，用户就可以对投影方向有更高的控制，以确保良好的网格质量。

 <img src="https://www.dealii.org/images/steps/developer/step-54.DirectionalProjection.png" alt="" width="500"> 


当然，后一种方法只有在表面的方向相当统一时才有效，这样就可以确定一个单一的投影方向。在表面方向接近投影方向的情况下，甚至有可能找不到方向性的投影。为了克服这些问题， OpenCASCADE::NormalToMeshProjectionBoundary 类实现了第三个投影算法。 OpenCASCADE::NormalToMeshProjectionBoundary 类的建立是将一个 <code>TopoDS_Shape</code> （至少包含一个面）分配给构造函数，其工作方式与 OpenCASCADE::DirectionalProjection. 完全一样。但是，正如该类的名字所暗示的， OpenCASCADE::NormalToMeshProjectionBoundary 试图想出一个合适的对要精化的网格元素的法线方向的估计，并将其用于新节点在CAD面上的投影。如果我们考虑二维空间中的网格边缘，其轴线方向是一个方向，沿着这个方向分割，以产生两个相同长度的新单元。我们在此将这一概念扩展到三维空间，并将所有新节点的投影方向近似于单元格的法线。

在下图中，受本教程中考虑的几何图形的启发，我们尝试比较所考虑的三种投影仪的行为。从左边可以看出，给定原始单元（蓝色），用法线投影找到的新点的位置不允许生成均匀的新元素（红色）。这种情况在进一步的细化步骤中会变得更糟。  由于我们考虑的几何体在某种程度上垂直于水平方向，以水平方向为投影方向的方向性投影（中心图像）在获得新的网格点方面做得相当好。然而，由于图片底部的表面几乎是水平的，我们可以预期在这些区域进行进一步细化步骤时，会出现问题。最后，右边的图片显示，位于单元轴上的节点将导致两个新单元具有相同的长度。当然，三维的情况会比这个简单的二维案例中描述的情况更复杂一些。然而，这个测试的结果证实，当考虑到任意形状的表面时，除非你有一个已知的更具体的方法，否则在测试的三种方法中，法线方向是最佳方法。


 <img src="https://www.dealii.org/images/steps/developer/step-54.ProjectionComparisons.png" alt="" width="700"> 




<a name="Thetestcase"></a><h3> The testcase </h3>


在这个程序中，我们将考虑为一个描述船头的真实几何体创建一个表面网格（这个几何体经常被用于CAD和网格生成的比较中，并且可以免费获得）。我们得到的表面网格可以用来解决边界元素方程，以模拟水在船舶周围的流动（类似于step-34的方式），但我们不会在这里尝试这样做。为了让你了解我们所考虑的几何形状，这里有一张图片。

 <img src="https://www.dealii.org/images/steps/developer/step-54.bare.png" alt="" width="500"> 

在程序中，我们从文件中读取几何体和粗略的网格，然后采用上面讨论的几个选项来放置新的顶点，进行一系列的网格细化步骤。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 
 * <a name="Includefiles"></a> 
 * <h3>Include files</h3>
 * 

 * 
 * We start with including a bunch of files that we will use in the
 * various parts of the program. Most of them have been discussed in
 * previous tutorials already:
 * 
 * @code
 * #include <deal.II/grid/tria.h>
 * #include <deal.II/grid/grid_generator.h>
 * #include <deal.II/grid/grid_in.h>
 * #include <deal.II/grid/grid_out.h>
 * #include <deal.II/numerics/data_out.h>
 * #include <deal.II/numerics/vector_tools.h>
 * 
 * @endcode
 * 
 * These are the headers of the opencascade support classes and
 * functions. Notice that these will contain sensible data only if you
 * compiled your deal.II library with support for OpenCASCADE, i.e.,
 * specifying <code>-DDEAL_II_WITH_OPENCASCADE=ON</code> and
 * <code>-DOPENCASCADE_DIR=/path/to/your/opencascade/installation</code>
 * when calling <code>cmake</code> during deal.II configuration.
 * 
 * @code
 * #include <deal.II/opencascade/manifold_lib.h>
 * #include <deal.II/opencascade/utilities.h>
 * 
 * 
 * @endcode
 * 
 * Finally, a few C++ standard header files
 * 
 * @code
 * #include <cmath>
 * #include <iostream>
 * #include <fstream>
 * #include <string>
 * 
 * @endcode
 * 
 * We isolate the rest of the program in its own namespace
 * 
 * @code
 * namespace Step54
 * {
 *   using namespace dealii;
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TheTriangulationOnCADclass"></a> 
 * <h3>The TriangulationOnCAD class</h3>
 * 

 * 
 * This is the main class. All it really does is store names for
 * input and output files, and a triangulation. It then provides
 * a function that generates such a triangulation from a coarse
 * mesh, using one of the strategies discussed in the introduction
 * and listed in the enumeration type at the top of the class.
 *   

 * 
 * The member functions of this class are similar to what you can
 * find in most of the other tutorial programs in the setup stage of
 * the grid for the simulations.
 * 

 * 
 * 
 * @code
 *   class TriangulationOnCAD
 *   {
 *   public:
 *     enum ProjectionType
 *     {
 *       NormalProjection       = 0,
 *       DirectionalProjection  = 1,
 *       NormalToMeshProjection = 2
 *     };
 * 
 * 
 *     TriangulationOnCAD(
 *       const std::string &  initial_mesh_filename,
 *       const std::string &  cad_file_name,
 *       const std::string &  output_filename,
 *       const ProjectionType surface_projection_kind = NormalProjection);
 * 
 *     void run();
 * 
 *   private:
 *     void read_domain();
 * 
 *     void refine_mesh();
 * 
 *     void output_results(const unsigned int cycle);
 * 
 *     Triangulation<2, 3> tria;
 * 
 *     const std::string initial_mesh_filename;
 *     const std::string cad_file_name;
 *     const std::string output_filename;
 * 
 *     const ProjectionType surface_projection_kind;
 *   };
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TriangulationOnCADTriangulationOnCAD"></a> 
 * <h4>TriangulationOnCAD::TriangulationOnCAD</h4>
 * 

 * 
 * The constructor of the TriangulationOnCAD class is very simple.
 * The input arguments are strings for the input and output file
 * names, and the enumeration type that determines which kind of
 * surface projector is used in the mesh refinement cycles (see
 * below for details).
 * 

 * 
 * 
 * @code
 *   TriangulationOnCAD::TriangulationOnCAD(
 *     const std::string &  initial_mesh_filename,
 *     const std::string &  cad_file_name,
 *     const std::string &  output_filename,
 *     const ProjectionType surface_projection_kind)
 *     : initial_mesh_filename(initial_mesh_filename)
 *     , cad_file_name(cad_file_name)
 *     , output_filename(output_filename)
 *     , surface_projection_kind(surface_projection_kind)
 *   {}
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TriangulationOnCADread_domain"></a> 
 * <h4>TriangulationOnCAD::read_domain</h4>
 * 

 * 
 * 

 * 
 * The following function represents the core of this program.  In
 * this function we import the CAD shape upon which we want to
 * generate and refine our triangulation. We assume that the CAD
 * surface is contained in the @p cad_file_name file (we provide an
 * example IGES file in the input directory called
 * "input/DTMB-5415_bulbous_bow.iges" that represents the bulbous bow of a
 * ship). The presence of several convex and concave high curvature
 * regions makes the geometry we provided a particularly meaningful
 * example.
 *   

 * 
 * After importing the hull bow surface, we extract some of the
 * curves and surfaces composing it, and use them to generate a set
 * of projectors. Such projectors define the rules the Triangulation
 * has to follow to position each new node during cell refinement.
 *   

 * 
 * To initialize the Triangulation, as done in previous tutorial
 * programs, we import a pre-existing grid saved in VTK format. We
 * assume here that the user has generated a coarse mesh
 * externally, which matches the IGES geometry. At the moment of
 * writing this tutorial, the
 * deal.II library does not automatically support generation of such
 * meshes, but there are several tools which can provide you with
 * reasonable initial meshes starting from CAD files.
 * In our example, the imported mesh is composed of a single
 * quadrilateral cell whose vertices have been placed on the CAD
 * shape.
 *   

 * 
 * After importing both the IGES geometry and the initial mesh, we
 * assign the projectors previously discussed to each of the edges
 * and cells which will have to be refined on the CAD surface.
 *   

 * 
 * In this tutorial, we will test the three different CAD surface
 * projectors described in the introduction, and will analyze the
 * results obtained with each of them.  As mentioned, each of these
 * projection strategies has been implemented in a different class,
 * and objects of these types can be assigned to a triangulation
 * using the Triangulation::set_manifold method.
 *   

 * 
 * The following function then first imports the given CAD file.
 * The function arguments are a string containing the desired file
 * name, and a scale factor. In this example, the scale factor is
 * set to 1e-3, as the original geometry is written in millimeters
 * (which is the typical unit of measure for most IGES files),
 * while we prefer to work in meters.  The output of the function
 * is an object of OpenCASCADE generic topological shape class,
 * namely a @p TopoDS_Shape.
 * 
 * @code
 *   void TriangulationOnCAD::read_domain()
 *   {
 *     TopoDS_Shape bow_surface = OpenCASCADE::read_IGES(cad_file_name, 1e-3);
 * 
 * @endcode
 * 
 * Each CAD geometrical object is defined along with a tolerance,
 * which indicates possible inaccuracy of its placement. For
 * instance, the tolerance @p tol of a vertex indicates that it can
 * be located in any point contained in a sphere centered in the
 * nominal position and having radius @p tol. While projecting a
 * point onto a surface (which will in turn have its tolerance) we
 * must keep in mind that the precision of the projection will be
 * limited by the tolerance with which the surface is built.
 * 

 * 
 * The following method extracts the tolerance of the given shape and
 * makes it a bit bigger to stay our of trouble:
 * 
 * @code
 *     const double tolerance = OpenCASCADE::get_shape_tolerance(bow_surface) * 5;
 * 
 * @endcode
 * 
 * We now want to extract a set of composite sub-shapes from the
 * generic shape. In particular, each face of the CAD file
 * is composed of a trimming curve of type @p TopoDS_Wire, which is
 * the collection of @p TopoDS_Edges that compose the boundary of a
 * surface, and a NURBS description of the surface itself. We will
 * use a line projector to associate the boundary of our
 * Triangulation to the wire delimiting the surface.  To extract
 * all compound sub-shapes, like wires, shells, or solids, we
 * resort to a method of the OpenCASCADE namespace.  The input of
 * OpenCASCADE::extract_compound_shapes is a shape and a set of empty
 * std::vectors of subshapes, which will be filled with all
 * compound shapes found in the given topological shape:
 * 
 * @code
 *     std::vector<TopoDS_Compound>  compounds;
 *     std::vector<TopoDS_CompSolid> compsolids;
 *     std::vector<TopoDS_Solid>     solids;
 *     std::vector<TopoDS_Shell>     shells;
 *     std::vector<TopoDS_Wire>      wires;
 * 
 *     OpenCASCADE::extract_compound_shapes(
 *       bow_surface, compounds, compsolids, solids, shells, wires);
 * 
 * @endcode
 * 
 * The next few steps are more familiar, and allow us to import an existing
 * mesh from an external VTK file, and convert it to a deal triangulation.
 * 
 * @code
 *     std::ifstream in;
 * 
 *     in.open(initial_mesh_filename);
 * 
 *     GridIn<2, 3> gi;
 *     gi.attach_triangulation(tria);
 *     gi.read_vtk(in);
 * 
 * @endcode
 * 
 * We output this initial mesh saving it as the refinement step 0.
 * 
 * @code
 *     output_results(0);
 * 
 * @endcode
 * 
 * The mesh imported has a single, two-dimensional cell located in
 * three-dimensional space. We now want to ensure that it is refined
 * according to the CAD geometry imported above. This this end, we get an
 * iterator to that cell and assign to it the manifold_id 1 (see
 * @ref GlossManifoldIndicator "this glossary entry").
 * We also get an iterator to its four faces, and assign each of them
 * the manifold_id 2:
 * 
 * @code
 *     Triangulation<2, 3>::active_cell_iterator cell = tria.begin_active();
 *     cell->set_manifold_id(1);
 * 
 *     for (const auto &face : cell->face_iterators())
 *       face->set_manifold_id(2);
 * 
 * @endcode
 * 
 * Once both the CAD geometry and the initial mesh have been
 * imported and digested, we use the CAD surfaces and curves to
 * define the projectors and assign them to the manifold ids just
 * specified.
 * 

 * 
 * A first projector is defined using the single wire contained in
 * our CAD file.  The ArclengthProjectionLineManifold will make
 * sure that every mesh edge located on the wire is refined with a
 * point that lies on the wire and splits it into two equal arcs
 * lying between the edge vertices. We first check
 * that the wires vector contains at least one element and then
 * create a Manifold object for it.
 *     

 * 
 * Once the projector is created, we then assign it to all the parts of
 * the triangulation with manifold_id = 2:
 * 
 * @code
 *     Assert(
 *       wires.size() > 0,
 *       ExcMessage(
 *         "I could not find any wire in the CAD file you gave me. Bailing out."));
 * 
 *     OpenCASCADE::ArclengthProjectionLineManifold<2, 3> line_projector(
 *       wires[0], tolerance);
 * 
 *     tria.set_manifold(2, line_projector);
 * 
 * @endcode
 * 
 * The surface projector is created according to what is specified
 * with the @p surface_projection_kind option of the constructor. In particular,
 * if the surface_projection_kind value equals @p NormalProjection, we select the
 * OpenCASCADE::NormalProjectionManifold. The new mesh points will
 * then initially be generated at the barycenter of the cell/edge
 * considered, and then projected on the CAD surface along its
 * normal direction.  The NormalProjectionManifold constructor
 * only needs a shape and a tolerance, and we then assign it to
 * the triangulation for use with all parts that manifold having id 1:
 * 
 * @code
 *     switch (surface_projection_kind)
 *       {
 *         case NormalProjection:
 *           {
 *             OpenCASCADE::NormalProjectionManifold<2, 3> normal_projector(
 *               bow_surface, tolerance);
 *             tria.set_manifold(1, normal_projector);
 * 
 *             break;
 *           }
 * 
 * @endcode
 * 
 * @p If surface_projection_kind value is @p DirectionalProjection, we select the
 * OpenCASCADE::DirectionalProjectionManifold class. The new mesh points
 * will then initially be generated at the barycenter of the cell/edge
 * considered, and then projected on the CAD surface along a
 * direction that is specified to the
 * OpenCASCADE::DirectionalProjectionManifold constructor. In this case,
 * the projection is done along the y-axis.
 * 
 * @code
 *         case DirectionalProjection:
 *           {
 *             OpenCASCADE::DirectionalProjectionManifold<2, 3>
 *               directional_projector(bow_surface,
 *                                     Point<3>(0.0, 1.0, 0.0),
 *                                     tolerance);
 *             tria.set_manifold(1, directional_projector);
 * 
 *             break;
 *           }
 * 
 * @endcode
 * 
 * As a third option, if @p surface_projection_kind value
 * is @p NormalToMeshProjection, we select the
 * OpenCASCADE::NormalToMeshProjectionManifold. The new mesh points will
 * again initially be generated at the barycenter of the cell/edge
 * considered, and then projected on the CAD surface along a
 * direction that is an estimate of the mesh normal direction.
 * The OpenCASCADE::NormalToMeshProjectionManifold constructor only
 * requires a shape (containing at least a face) and a
 * tolerance.
 * 
 * @code
 *         case NormalToMeshProjection:
 *           {
 *             OpenCASCADE::NormalToMeshProjectionManifold<2, 3>
 *               normal_to_mesh_projector(bow_surface, tolerance);
 *             tria.set_manifold(1, normal_to_mesh_projector);
 * 
 *             break;
 *           }
 * 
 * @endcode
 * 
 * Finally, we use good software cleanliness by ensuring that this
 * really covers all possible options of the @p case statement. If we
 * get any other value, we simply abort the program:
 * 
 * @code
 *         default:
 *           AssertThrow(false, ExcInternalError());
 *       }
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TriangulationOnCADrefine_mesh"></a> 
 * <h4>TriangulationOnCAD::refine_mesh</h4>
 * 

 * 
 * This function globally refines the mesh. In other tutorials, it
 * would typically also distribute degrees of freedom, and resize
 * matrices and vectors. These tasks are not carried out here, since
 * we are not running any simulation on the Triangulation produced.
 *   

 * 
 * While the function looks innocent, this is where most of the work we are
 * interested in for this tutorial program actually happens. In particular,
 * when refining the quads and lines that define the surface of the ship's
 * hull, the Triangulation class will ask the various objects we have
 * assigned to handle individual manifold ids for where the new vertices
 * should lie.
 * 
 * @code
 *   void TriangulationOnCAD::refine_mesh()
 *   {
 *     tria.refine_global(1);
 *   }
 * 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TriangulationOnCADoutput_results"></a> 
 * <h4>TriangulationOnCAD::output_results</h4>
 * 

 * 
 * Outputting the results of our computations is a rather mechanical
 * task. All the components of this function have been discussed
 * before:
 * 
 * @code
 *   void TriangulationOnCAD::output_results(const unsigned int cycle)
 *   {
 *     const std::string filename =
 *       (output_filename + "_" + Utilities::int_to_string(cycle) + ".vtk");
 *     std::ofstream logfile(filename);
 *     GridOut       grid_out;
 *     grid_out.write_vtk(tria, logfile);
 *   }
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="TriangulationOnCADrun"></a> 
 * <h4>TriangulationOnCAD::run</h4>
 * 

 * 
 * This is the main function. It should be self explanatory in its
 * briefness:
 * 
 * @code
 *   void TriangulationOnCAD::run()
 *   {
 *     read_domain();
 * 
 *     const unsigned int n_cycles = 5;
 *     for (unsigned int cycle = 0; cycle < n_cycles; ++cycle)
 *       {
 *         refine_mesh();
 *         output_results(cycle + 1);
 *       }
 *   }
 * } // namespace Step54
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h3>The main() function</h3>
 * 

 * 
 * This is the main function of this program. It is in its basic structure
 * like all previous tutorial programs, but runs the main class through the
 * three possibilities of new vertex placement:
 * 
 * @code
 * int main()
 * {
 *   try
 *     {
 *       using namespace Step54;
 * 
 *       const std::string in_mesh_filename = "input/initial_mesh_3d.vtk";
 *       const std::string cad_file_name    = "input/DTMB-5415_bulbous_bow.iges";
 * 
 *       std::cout << "----------------------------------------------------------"
 *                 << std::endl;
 *       std::cout << "Testing projection in direction normal to CAD surface"
 *                 << std::endl;
 *       std::cout << "----------------------------------------------------------"
 *                 << std::endl;
 *       std::string        out_mesh_filename = ("3d_mesh_normal_projection");
 *       TriangulationOnCAD tria_on_cad_norm(in_mesh_filename,
 *                                           cad_file_name,
 *                                           out_mesh_filename,
 *                                           TriangulationOnCAD::NormalProjection);
 *       tria_on_cad_norm.run();
 *       std::cout << "----------------------------------------------------------"
 *                 << std::endl;
 *       std::cout << std::endl;
 *       std::cout << std::endl;
 * 
 *       std::cout << "----------------------------------------------------------"
 *                 << std::endl;
 *       std::cout << "Testing projection in y-axis direction" << std::endl;
 *       std::cout << "----------------------------------------------------------"
 *                 << std::endl;
 *       out_mesh_filename = ("3d_mesh_directional_projection");
 *       TriangulationOnCAD tria_on_cad_dir(
 *         in_mesh_filename,
 *         cad_file_name,
 *         out_mesh_filename,
 *         TriangulationOnCAD::DirectionalProjection);
 *       tria_on_cad_dir.run();
 *       std::cout << "----------------------------------------------------------"
 *                 << std::endl;
 *       std::cout << std::endl;
 *       std::cout << std::endl;
 * 
 *       std::cout << "----------------------------------------------------------"
 *                 << std::endl;
 *       std::cout << "Testing projection in direction normal to mesh elements"
 *                 << std::endl;
 *       std::cout << "----------------------------------------------------------"
 *                 << std::endl;
 *       out_mesh_filename = ("3d_mesh_normal_to_mesh_projection");
 *       TriangulationOnCAD tria_on_cad_norm_to_mesh(
 *         in_mesh_filename,
 *         cad_file_name,
 *         out_mesh_filename,
 *         TriangulationOnCAD::NormalToMeshProjection);
 *       tria_on_cad_norm_to_mesh.run();
 *       std::cout << "----------------------------------------------------------"
 *                 << std::endl;
 *       std::cout << std::endl;
 *       std::cout << std::endl;
 *     }
 *   catch (std::exception &exc)
 *     {
 *       std::cerr << std::endl
 *                 << std::endl
 *                 << "----------------------------------------------------"
 *                 << std::endl;
 *       std::cerr << "Exception on processing: " << std::endl
 *                 << exc.what() << std::endl
 *                 << "Aborting!" << std::endl
 *                 << "----------------------------------------------------"
 *                 << std::endl;
 * 
 *       return 1;
 *     }
 *   catch (...)
 *     {
 *       std::cerr << std::endl
 *                 << std::endl
 *                 << "----------------------------------------------------"
 *                 << std::endl;
 *       std::cerr << "Unknown exception!" << std::endl
 *                 << "Aborting!" << std::endl
 *                 << "----------------------------------------------------"
 *                 << std::endl;
 *       return 1;
 *     }
 * 
 *   return 0;
 * }
 * @endcode
examples/step-54/doc/results.dox



<a name="Results"></a><h1>Results</h1>


程序的执行会产生一系列的网格文件 <code>3d_mesh_*.vtk</code> ，我们可以用任何可以读取VTK文件格式的常用可视化程序来进行可视化。

下表说明了采用正常投影策略得到的结果。表中的前两行显示的是逐步细化的网格的侧视图，覆盖在精确几何体的非常精细的渲染上。深红色和浅红色的区域只是表示当前的网格或精细的几何体更接近观察者；这种区别没有任何特别深刻的意义。最后一排图片描述了第二排中相同网格的正视图（镜像到几何体的两边）。


 <table style="width:90%">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.common_0.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_1.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_2.png" alt="" width="400"></td>
  </tr>
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_3.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_4.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_5.png" alt="" width="400"></td>
  </tr>
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_front_3.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_front_4.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_front_5.png" alt="" width="400"></td>
  </tr>
</table> 

从图片中可以看出，正如我们所预料的那样，当应用于具有明显曲率变化的表面时，正常的细化策略无法产生良好的形状的元素。这在船体的球体上尤其明显，所有的新点都被放置在球体的上部，而下部则完全没有被解决。

下表的排列方式与上表相同，说明了采用方向性投影方法获得的结果，其中选择的投影方向是Y轴（在每幅图像的左下方用一个小的黄色箭头表示）。


 <table style="width:90%">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.common_0.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.directional_1.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.directional_2.png" alt="" width="400"></td>
  </tr>
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.directional_3.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.directional_4.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.directional_5.png" alt="" width="400"></td>
  </tr>
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.directional_front_3.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.directional_front_4.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.directional_front_5.png" alt="" width="400"></td>
  </tr>
</table> 

这些图像证实，用定向投影得到的网格质量明显高于沿表面法线投影得到的网格。然而，在球体底部观察到一些在Y方向上拉长的元素，那里的表面几乎与选择的投影方向平行。

最后的测试显示了使用面的法线投影的结果。

 <table style="width:90%">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.common_0.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_to_mesh_1.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_to_mesh_2.png" alt="" width="400"></td>
  </tr>
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_to_mesh_3.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_to_mesh_4.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_to_mesh_5.png" alt="" width="400"></td>
  </tr>
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_to_mesh_front_3.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_to_mesh_front_4.png" alt="" width="400"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-54.normal_to_mesh_front_5.png" alt="" width="400"></td>
  </tr>
</table> 

图片证实了法线投影的方法导致网格在整个细化步骤中保持均匀的间隔。同时，这些网格很好地表现了原始的几何形状，甚至在灯泡的底部区域也是如此，这一点在使用定向投影仪或法线投影仪时并没有得到很好的恢复。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-54.cc"
*/
