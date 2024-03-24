#include "Mesh.H"

/* This file defines an unstructured mesh class to be used for CFD purposes.
At the base level, the mesh is defined by a set of points (defined in supplied
pointsFile). Faces are made up of specific sets of points as defined in the
facesFile. Volumetric cells are made up of a set of faces as defined in
cellsFile. Boundary patches are made up of closed sets of faces (as defined in
boundaryFile).

The Mesh class communicates with four other classes: Point, Face, Cell, and
Boundarypatch. A mesh has total control over its vectors of Points, Faces,
Cells, and Boundarypatches in that these vectors are private to the mesh and the
mesh may modify the instances in these vectors.

All grids generated will be inherently 3D. We can obtain equivalents
for 1D and 2D models by restricting the number of cells in the z-direction (2D)
and y-direction (1D) to 1. Faces pointing in these directions will be designated
as "empty" (no boundary condition necessary).

Note that a significant amount of parsing is necessary when reading the files,
as the format is string-based (.dat files)
*/

// Constructor to read mesh data from files
Mesh::Mesh(const std::string& pointsFile, const std::string& facesFile,
           const std::string& cellsFile, const std::string& boundaryFile) {
  // Read files containing specifications of mesh
  readPoints(pointsFile);
  readFaces(facesFile);
  readCells(cellsFile);
  readBoundaryPatches(boundaryFile);

  // Calculate area-scaled normal vectors to cell faces
  // and face centroids
  // Assign these to the relevant face instances
  for (int i = 0; i < faces.size(); i++) {
    Set_Sf(i);
    set_face_centroid(i);
  }

  // Calculate cell centroids and cell volumes
  // Assign these to the relevant cell instances
  for (int i = 0; i < cells.size(); i++) {
    set_cell_centroid(i);
    set_cell_vol(i);
  }

  // Assign face owner and neighbour cell (i.e. the cells of whose boundary the
  // face is part). We distinguish owner and neighbour cells for
  // consistency purposes when solving associated system of equations. Face
  // owner is the cell for which the face normal is outward-pointing
  assign_face_owners();

  // For each cell assign a set of neighbouring cells
  assign_cell_neighbours();

  /* For each boundary patch assign to it an owner cell
  // This function also ensures that face normals are outward
   pointing for faces that form part of a boundary patch */
  assign_boundary_cells();

  /* Face interpolation factor and delta coefficient
  are attributes of each face. Only assign face-interpolation
  factor for non-boundary faces */
  for (int i = 0; i < faces.size(); i++) {
    if (faces[i].boundary < 0) {
      set_face_interpolation_factor(i);
    }
    set_delta_coef(i);
  }
}

/* Function to read points from a file. The pointsFile will always have the
form: nrPoints (number of points)
(
(p1_x p1_y p1_z)
(p2_x p2_y p2_z)
...
(pn_x pn_y pn_z)
)
*/
void Mesh::readPoints(const std::string& fileName) {
  std::ifstream inputFile(fileName);

  if (!inputFile.is_open()) {
    std::cerr << "Failed to open the points file." << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string line;     // declare string to hold current line
  std::string str_num;  // declare string to representation of current number
  std::getline(inputFile, line);  // get first line from file
  int numPoints;                  // will hold number of vertices

  // Get string representation of nr of points
  for (std::size_t i = 0; i < line.length(); ++i) {
    str_num += line[i];
  }

  // convert to integer
  numPoints = std::stoi(str_num);

  // Resize points vector (vector of points objects declared in header file)
  points.resize(numPoints);

  std::getline(inputFile, line);  // skip past line containing only parenthesis

  int j;  // declare variable that will be used to iterate over characters in
          // line

  for (int i = 0; i < numPoints;
       ++i) {  // iterate over each line containing vertex
    std::getline(inputFile, line);

    j = 1;  // Each line starts with a parenthesis that we can skip past
    str_num = "";
    while (line[j] != ' ') {
      str_num += line[j];  // append characters of number to string until space
                           // is reached
      j += 1;
    }

    points[i].coord[0] =
        std::stod(str_num);  // Assign to given point its x-coordinate

    j += 1;  // skip past space character
    str_num = "";
    while (line[j] != ' ') {
      str_num += line[j];
      j += 1;
    }

    points[i].coord[1] =
        std::stod(str_num);  // Assign to given point its y-coordinate

    j += 1;
    str_num = "";
    while ((line[j] != ' ') && (line[j] != ')')) {
      str_num += line[j];
      j += 1;
    }

    points[i].coord[2] =
        std::stod(str_num);  // Assign to given point its z-coordinate
  }

  inputFile.close();
}

/* Function to read face definitions from a file.
Let k_i denote the number of points that make up the corners of a face
The facesFile will always have the form of:
nrFaces (number of faces)
(
k_0([space separated list of k_0 points (indexed from zero) for face 0])
k_1([space separated list of k_1 points (indexed from zero) for face 1])
...
k_n ([space separated list of k_n points (indexed from zero) for face n])
)

**Note: it is assumed that points are provided in a counterclockwise order
with respect to the face plane
*/
void Mesh::readFaces(const std::string& fileName) {
  std::ifstream inputFile(fileName);

  if (!inputFile.is_open()) {
    std::cerr << "Failed to open the faces file." << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string line;     // declare string to hold current line
  std::string str_num;  // declare string to representation of current number
  std::getline(inputFile, line);  // get first line from file
  int numFaces;                   // will hold number of faces
  int numVertices;  // will hold number of verteces on current face

  // Get string representation of nr of faces
  for (std::size_t i = 0; i < line.length(); ++i) {
    str_num += line[i];
  }

  // convert string to integer and resize face vector of face objects (decalred
  // in header)
  numFaces = std::stoi(str_num);
  faces.resize(numFaces);

  std::getline(inputFile, line);  // skip past line containing only parenthesis

  int j;  // declare variable that will be used to iterate over characters in
          // line

  /*Iteration over each line proceeds much as with the points, but now we
  populate the "vertices"-vector for each face object */
  for (int i = 0; i < numFaces;
       ++i) {  // iterate over each line containing list of vertices
    j = 0;     // reset j
    std::getline(inputFile, line);
    str_num = "";

    // Get number of points which make up a face in string form
    while ((line[j] != '(') && (line[j] != ' ')) {
      str_num += line[j];
      j += 1;
    }

    // Convert to integer and resize vertices vector of the given face object
    numVertices = std::stoi(str_num);
    faces[i].vertices.resize(numVertices);

    j += 1;  // skip to first character of first vertex label on this face

    for (int k = 0; k < numVertices; ++k) {
      str_num = "";
      while ((line[j] != ' ') && (line[j] != ')')) {
        str_num += line[j];
        j += 1;
      }

      faces[i].vertices[k] = std::stoi(str_num);
      j += 1;
    }

    /* Faces contain integer attribute boundary which indicates the index of the
    patch if positive and is -1 if the face is not on boundary of domain

    They also contain boolean attribute empty_face which indicate if face
    is on a boundary that is being considered with the boundary conditions or
    not (only relevant for 1D and 2D models)

    Initially mark all faces as not being on a boundary and not being empty
    This will be changed when boundary patches are initialised and
    when cell-neighbours are assigned, respectively.
    A face is empty if it is at the boundary of mesh but forms part
    of no boundary patch.
    */
    faces[i].boundary = -1;
    faces[i].empty_face = false;
  }

  inputFile.close();
}

/* Function to read cells from a file
Let k_i denote the number of faces that make up the boundaries of a cell
The cellsFile will always have the form of:
nrcells (number of cells)
(
k_0([space separated list of k_0 faces (indexed from zero) for cell 0])
k_1([space separated list of k_1 faces (indexed from zero) for cell 1])
...
k_n ([space separated list of k_n faces (indexed from zero) for cell n])
)
*/
void Mesh::readCells(const std::string& fileName) {
  std::ifstream inputFile(fileName);

  if (!inputFile.is_open()) {
    std::cerr << "Failed to open the cells file." << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string line;     // declare string to hold current line
  std::string str_num;  // declare string to representation of current number
  std::getline(inputFile, line);  // get first line from file
  int numCells;                   // will hold number of cells
  int numFaces;  // will hold number of verteces on current Cell

  // get number of cells
  for (std::size_t i = 0; i < line.length(); ++i) {
    str_num += line[i];
  }

  // Resize vector that holds cell objects (declared in header file)
  numCells = std::stoi(str_num);
  cells.resize(numCells);

  std::getline(inputFile, line);  // skip past line containing only parenthesis

  int j;  // declare variable that will be used to iterate over characters in
          // line

  /* Note that we cannot directly assign to faces their owner and neighbour
  cells at this stage, as we do not know what faces belong to boundary patches.
  Thus we assign to every face an "owner" and "neighbour" value of -1 */

  for (int i = 0; i < numCells;
       ++i) {  // iterate over each line containing list of faces for each cell
    j = 0;     // reset j
    std::getline(inputFile, line);
    str_num = "";

    // find number of faces belonging to current cell
    while ((line[j] != '(') && (line[j] != ' ')) {
      str_num += line[j];
      j += 1;
    }

    // Resize faces vector of current cell object
    numFaces = std::stoi(str_num);
    cells[i].faces.resize(numFaces);

    j += 1;  // skip to first character of first face index on this cell

    // Iterate over number of faces defining current cell
    for (int k = 0; k < numFaces; ++k) {
      str_num = "";
      while ((line[j] != ' ') && (line[j] != ')')) {
        str_num += line[j];
        j += 1;
      }

      // Populate faces vector in current cell object with current face index
      // initially define every neighbour and owner by -1
      // which will signify boundary patch
      cells[i].faces[k] = std::stoi(str_num);
      faces[cells[i].faces[k]].ownerCell = -1;
      faces[cells[i].faces[k]].neighbourCell = -1;
      j += 1;
    }
  }

  inputFile.close();
}

/* Function to read boundary patches from a file
Let k_i denote the number of faces that make up a boundary patch
Let name_i denote the name of patch i
The boundaryFile will always have the form of:
nrpatches (number of boundary patches)
(
name_0
k_0
(
[space separated list of k_0 faces (indexed from zero) for patch 0])
)
...
name_n
k_n
(
([space separated list of k_n faces (indexed from zero) for patch n])
)
)
*/
void Mesh::readBoundaryPatches(const std::string& fileName) {
  std::ifstream inputFile(fileName);

  if (!inputFile.is_open()) {
    std::cerr << "Failed to open the boundary patches file." << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string line;     // declare string to hold current line
  std::string str_num;  // declare string to representation of current number
  int numBoundaries;    // will hold number of Boundary patches
  int numFaces;  // will hold number of verteces on current Boundary patch

  std::getline(inputFile, line);  // get first line from file

  // get number of boundary patches
  for (std::size_t i = 0; i < line.length(); ++i) {
    str_num += line[i];
  }

  numBoundaries = std::stoi(str_num);
  boundaryPatches.resize(numBoundaries);

  std::getline(inputFile, line);  // jump to line containing only parenthesis

  int j;  // declare variable that will be used to iterate over characters in
          // line

  for (int i = 0; i < numBoundaries;
       i++) {  // iterate over each block containing name of patch and info on
               // patches
    j = 0;     // reset j
    std::getline(inputFile, line);  // jump to line containing name

    boundaryPatches[i].name = line;  // assign name

    std::getline(inputFile,
                 line);  // jump to line containing number of faces adjoining

    str_num = "";

    // find number of faces adjoining current boundary patch
    while (line[j]) {
      str_num += line[j];
      j += 1;
    }

    numFaces = std::stoi(str_num);
    boundaryPatches[i].faces.resize(numFaces);

    std::getline(inputFile,
                 line);  // jump to line containing which faces are adjoining
    std::getline(inputFile, line);
    j = 0;  // Reset j

    // for each boundary patch, assign to its faces vector the right face labels
    for (int k = 0; k < numFaces; ++k) {
      str_num = "";
      while ((line[j] != ' ') && (line[j] != ')') && (j < line.length())) {
        str_num += line[j];
        j += 1;
      }
      boundaryPatches[i].faces[k] = std::stoi(str_num);

      // Assign to each face adjoining patch the patch number
      faces[std::stoi(str_num)].boundary = i;

      j += 1;
    }
    std::getline(inputFile, line);
  }

  inputFile.close();
}

void Mesh::assign_face_owners() {
  // Need to loop over cells to assign face owners and neighbours
  // as faces do not point to cells in any way
  int numCells = cells.size();
  Eigen::ArrayXi cell_faces;
  Eigen::Vector3d f_centroid, Sf, c_centroid;
  /*Declare vectors that will hold face areal centroid,
  area-scaled normal, and cell-centroid.  */

  // Loop over cells
  for (int i = 0; i < numCells; i++) {
    // list of faces for each cell
    Eigen::ArrayXi cell_faces = cells[i].faces;
    int numFaces_cell = cell_faces.size();

    // centroid of cell
    c_centroid = cells[i].centroid;

    /*For each face adjoining a cell, if dot product of
    (face centroid - cell centroid) with area-scaled normal
     is positive, cell is owner of the face */
    for (int j = 0; j < numFaces_cell; j++) {
      // Get area-scaled face normal
      Eigen::Vector3d Sf = faces[cell_faces[j]].Sf;

      // Get face centroid
      f_centroid = faces[cell_faces[j]].centroid;

      // If face is internal
      if (faces[cell_faces[j]].boundary < 0) {
        double dot = Sf.dot(f_centroid - c_centroid);

        if (dot > 0) {
          faces[cell_faces[j]].ownerCell = i;
        }

        else {
          faces[cell_faces[j]].neighbourCell = i;
        }
      }

      // if face is at boundary, adjoining cell is automatically owner
      else {
        faces[cell_faces[j]].ownerCell = i;
      }
    }
  }

  // Now make sure each face has an owner. That is, if a face only has a
  // neighbour cell but no owner, make the neighbour the owner and set neighbour
  // to -1
  for (int j = 0; j < faces.size(); j++) {
    if (faces[j].ownerCell < 0) {
      faces[j].ownerCell = faces[j].neighbourCell;
      faces[j].neighbourCell = -1;
    }
  }
}

void Mesh::assign_cell_neighbours() {
  // loop over faces and couple face owner/neighbours as cell neighbours
  for (int i = 0; i < faces.size(); i++) {
    int cell1, cell2;
    cell1 = faces[i].ownerCell;
    cell2 = faces[i].neighbourCell;

    if ((cell2 < 0) && (faces[i].boundary < 0)) {
      // if face is on boundary and is not already assigned to a patch
      //(done during parsing), mark it as empty face.
      faces[i].empty_face = true;
    }

    // Else face marks boundary between two cells. Alter the cells'
    // neighbour vectors accordingly
    else if (faces[i].boundary < 0) {
      cells[cell1].nbrs.push_back(cell2);
      cells[cell2].nbrs.push_back(cell1);
    }
  }
}

void Mesh::assign_boundary_cells() {
  // Loop over boundary patches and assign adjoining cell as owner of first face
  int a_face, face_nr, cell_nr, numFaces;
  Eigen::Vector3d cell_cent, face_cent, Sf, f_normal;
  double distanceSq_add_normal, distanceSq_sub_normal;

  for (int i = 0; i < boundaryPatches.size(); i++) {
    a_face = boundaryPatches[i].faces[0];
    boundaryPatches[i].cell = faces[a_face].ownerCell;

    // also make sure that Sf (area-scaled normal) is outward pointing for every
    // face adjoining a boundary patch:
    numFaces = boundaryPatches[i].faces.size();

    for (int j = 0; j < numFaces; j++) {
      /* Procedure is the same as when owner/neighbour cells of faces
       were assigned. If Sf is inward pointing, we have
       (face centroid - cell centroid) dot Sf < 0*/
      face_nr = boundaryPatches[i].faces[j];
      cell_nr = faces[face_nr].ownerCell;
      cell_cent = cells[cell_nr].centroid;
      face_cent = faces[face_nr].centroid;
      Sf = faces[face_nr].Sf;

      if ((face_cent - cell_cent).dot(Sf) < 0) {
        faces[face_nr].Sf *= -1;
      }
    }
  }
}

// Add getter and setter functions for various mesh properties
int Mesh::getNumPoints() const { return points.size(); }

int Mesh::getNumFaces() const { return faces.size(); }

int Mesh::getNumCells() const { return cells.size(); }

int Mesh::getNumBoundaryPatches() const { return boundaryPatches.size(); }

std::vector<int> Mesh::get_cell_nbrs(int cell_nr) const {
  return cells[cell_nr].nbrs;
}

double Mesh::get_cell_vol(int cell_nr) const { return cells[cell_nr].vol; }

Eigen::Vector3d Mesh::get_cell_velocity(int cell_nr) const {
  return cells[cell_nr].vel;
}

bool Mesh::get_face_empty(int face_nr) const {
  return faces[face_nr].empty_face;
}

int Mesh::get_face_owner(int face_nr) const { return faces[face_nr].ownerCell; }

int Mesh::get_face_nbr(int face_nr) const {
  return faces[face_nr].neighbourCell;
}

Eigen::Vector3d Mesh::get_face_normal(int face_nr) const {
  return faces[face_nr].Sf;
}

int Mesh::face_boundary(int face_nr) const { return faces[face_nr].boundary; }

Eigen::Vector3d Mesh::get_cell_centroid(int cell_nr) const {
  return cells[cell_nr].centroid;
}

double Mesh::get_face_interpolation_factor(int face_nr) const {
  return faces[face_nr].f_x;
}

double Mesh::get_delta_coef(int face_nr) const { return faces[face_nr].delta; }

// Get the gradient of a field at a cell
const Eigen::Vector3d Mesh::get_cell_grad(int cell_nr, Field field) const {
  switch (field) {
    case Field::vx:
      return cells[cell_nr].grad_vx;
    case Field::vy:
      return cells[cell_nr].grad_vy;
    case Field::vz:
      return cells[cell_nr].grad_vz;
    case Field::p:
      return cells[cell_nr].grad_p;
  }
}

// Set gradient in a cell
void Mesh::set_cell_grad(int cell_nr, Field field, Eigen::Vector3d grad) {
  switch (field) {
    case Field::vx:
      cells[cell_nr].grad_vx = grad;
      return;
    case Field::vy:
      cells[cell_nr].grad_vy = grad;
      return;
    case Field::vz:
      cells[cell_nr].grad_vz = grad;
      return;
    case Field::p:
      cells[cell_nr].grad_p = grad;
      return;
  }
}

// Get the provided field value at cell centroid of the given cell
double Mesh::get_cell_val(int cell_nr, Field field) const {
  switch (field) {
    case Field::vx:
      return cells[cell_nr].vel(0);
    case Field::vy:
      return cells[cell_nr].vel(1);
    case Field::vz:
      return cells[cell_nr].vel(2);
    case Field::p:
      return cells[cell_nr].p;
  }
}

// Get the net flux out of a cell. This is defined as the flux
// over each of the cell's faces. If cell is neighbour to a face
// the sign of the flux over that face must be inverted
double Mesh::get_cell_netflux(int cell_nr) const {
  double cell_flux = 0.0;
  for (auto idx : cells[cell_nr].faces) {
    if (!get_face_empty(idx)) {
      int owner = faces[idx].ownerCell;

      if (cell_nr == owner) {
        cell_flux += faces[idx].flux;
      }

      else {
        cell_flux -= faces[idx].flux;
      }
    }
  }
  return cell_flux;
}

// Return an Eigen array containing the given field values at each cell
const Eigen::ArrayXd Mesh::get_field(Field field) const {
  Eigen::ArrayXd field_array(cells.size());

  switch (field) {
    case Field::vx:
      for (int i = 0; i < cells.size(); i++) {
        field_array(i) = cells[i].vel(0);
      }
      return field_array;
    case Field::vy:
      for (int i = 0; i < cells.size(); i++) {
        field_array(i) = cells[i].vel(1);
      }
      return field_array;
    case Field::vz:
      for (int i = 0; i < cells.size(); i++) {
        field_array(i) = cells[i].vel(1);
      }
      return field_array;
    case Field::p:
      for (int i = 0; i < cells.size(); i++) {
        field_array(i) = cells[i].p;
      }
      return field_array;
  }
}

void Mesh::set_velocity(Eigen::MatrixXd& v) {
  for (int i = 0; i < cells.size(); i++) {
    set_cell_velocity(i, v(i, 0), v(i, 1), v(i, 2));
  }
}

void Mesh::set_cell_velocity(int cell_nr, double vx, double vy, double vz) {
  cells[cell_nr].vel(0) = vx;
  cells[cell_nr].vel(1) = vy;
  cells[cell_nr].vel(2) = vz;
}

void Mesh::set_pressure(Eigen::ArrayXd& p) {
  for (int i = 0; i < cells.size(); i++) {
    set_cell_pressure(i, p(i));
  }
}

void Mesh::set_cell_pressure(int cell_nr, double p) { cells[cell_nr].p = p; }

// Return volume of all cells
Eigen::VectorXd Mesh::getCellVolumes() const {
  int numCells = cells.size();
  Eigen::VectorXd volumes(numCells);

  for (int i = 0; i < numCells; i++) {
    volumes[i] = cells[i].vol;
  }

  return volumes;
}

// face_raw_avg is used to find the raw average of the points
// defining face corners
Eigen::Vector3d Mesh::face_raw_avg(int face_nr) const {
  // find number of vertices for current face
  int numVertices = faces[face_nr].vertices.size();

  // declare vector that will hold RAW AVERAGE of all vertices
  // and integer variable that will hold a point index
  Eigen::Vector3d raw_avg(0.0, 0.0, 0.0);
  int vertex_lable;

  // compute raw average
  for (int i = 0; i < numVertices; i += 1) {
    vertex_lable = faces[face_nr].vertices[i];
    raw_avg += (points[vertex_lable].coord / numVertices);
  }

  return raw_avg;
}

/*Set Sf is used to compute face normal scaled by face area. Face will not be a
perfect plane due to numerical inaccuracy. Hence, calculate face area by summing
areas of triangles formed by adjacent points and face raw average (on a given
face). The DIRECTION of Sf is found by summing the normal vectors of said
triangles. This sum is then normalised and subsequently scaled by the total
area*/
void Mesh::Set_Sf(int face_nr) {
  std::vector<std::array<int, 2> > tri_list = triangle_list(face_nr);
  Eigen::Vector3d raw_avg = face_raw_avg(face_nr);
  Eigen::Vector3d Sf(0, 0, 0);

  // Variable to contain total face area
  double areas = 0;

  // Define two "vectors" that will be used to compute normal to
  // each triangle of the face
  // In practice just subtract third vertex from first and second
  Eigen::Vector3d vec1, vec2;

  /*tri_list contains sets of pairs of adjacent points/face corners
  two points will only form a pair if they are adjacent when going
  round the face counter-clockwise*/
  for (int i = 0; i < tri_list.size(); i++) {
    // Vectors from each of points in a pair to raw face average
    vec1 = points[tri_list[i][0]].coord - raw_avg;
    vec2 = points[tri_list[i][1]].coord - raw_avg;

    // Vector normal to triangle formed by pair and face raw avg
    Sf += vec1.cross(vec2);

    // Add to total area of face
    areas += vec1.cross(vec2).norm() / 2;
  }

  Sf.normalize();
  Sf = areas * Sf;

  // Set Sf attribute of given face
  faces[face_nr].Sf = Sf;
}

/*triangle_list produces a vector of pairs of adjacent points/face corners.
Two points will only form a pair if they are adjacent when going
round the face counter-clockwise*/
std::vector<std::array<int, 2> > Mesh::triangle_list(int face_nr) const {
  Eigen::Vector3d raw_avg = face_raw_avg(face_nr);

  // find number of vertices for current face
  int numVertices = faces[face_nr].vertices.size();

  // declare array that will hold list of vertices for non-overlapping triangles
  // for face where vertex-points are identified by their indices. Only two
  // vertices are listed for each triangle since the third will always be the
  // raw_avg
  std::vector<std::array<int, 2> > tri_list;
  tri_list.resize(numVertices);

  for (int i = 0; i < numVertices - 1; i++) {
    tri_list[i][0] = faces[face_nr].vertices[i];
    tri_list[i][1] = faces[face_nr].vertices[i + 1];
  }

  // Last vertex must be coupled with first
  tri_list[numVertices - 1][0] = faces[face_nr].vertices[numVertices - 1];
  tri_list[numVertices - 1][1] = faces[face_nr].vertices[0];

  return tri_list;
}

// Find midpoint of each triangle of which face is composed.
// This function will be used when calculating face centroid
std::vector<Eigen::Vector3d> Mesh::triangle_midpoints(int face_nr) const {
  // Last vertex is always the raw average of all vertices
  std::vector<std::array<int, 2> > tri_list = triangle_list(face_nr);
  int numVertices = faces[face_nr].vertices.size();
  Eigen::Vector3d raw_avg = face_raw_avg(face_nr);
  std::vector<Eigen::Vector3d> mpoints;
  mpoints.resize(numVertices);

  // midpoint of each triangle is just average of vertices
  for (int i = 0; i < numVertices; i++) {
    mpoints[i] = (points[tri_list[i][0]].coord + points[tri_list[i][1]].coord +
                  raw_avg) /
                 3.0;
  }

  return mpoints;
}

void Mesh::set_face_centroid(int face_nr) {
  // Compute face centroid by taking weighted mean of midpoints
  Eigen::Vector3d centroid(0.0, 0.0, 0.0);
  std::vector<std::array<int, 2> > tri_list = triangle_list(face_nr);
  int numVertices = faces[face_nr].vertices.size();
  std::vector<Eigen::Vector3d> midpoints = triangle_midpoints(face_nr);
  Eigen::Vector3d raw_avg = face_raw_avg(face_nr);

  // Define two "vectors" that will be used to compute normal to
  // each triangle of the face
  // In practice just subtract third vertex from first and second
  Eigen::Vector3d vec1, vec2;

  // Define variable for area of given triangle in face and variable for
  // total area
  double tri_area, total_area = 0;

  // Loop over pairs in tri_list and calculated midpoints to compute
  // weighted average
  for (int i = 0; i < numVertices; i++) {
    // Vectors from each of points in a pair to raw face average
    vec1 = points[tri_list[i][0]].coord - raw_avg;
    vec2 = points[tri_list[i][1]].coord - raw_avg;

    // area of triangle
    tri_area = vec1.cross(vec2).norm();
    total_area += tri_area;

    centroid += (tri_area * midpoints[i]);
  }

  centroid /= total_area;

  faces[face_nr].centroid = centroid;
}

// Find raw average of face centroids around cell
// This is just a simple average where we loop over faces
// adjoining cell
Eigen::Vector3d Mesh::cell_raw_avg(int cell_nr) const {
  int numFaces = cells[cell_nr].faces.size();
  Eigen::ArrayXi cell_faces = cells[cell_nr].faces;
  Eigen::Vector3d raw_avg(0.0, 0.0, 0.0), face_cent;

  for (int i = 0; i < numFaces; i++) {
    face_cent = faces[cell_faces[i]].centroid;
    raw_avg += face_cent;
  }

  raw_avg /= numFaces;

  return raw_avg;
}

/* Cell centroid will be calculated as the average of face-centroids,
weighted by the volume of the pyramid formed between face and cell raw average.
The given function returns a vector of pyramid volumes for a given cell, indexed
in order of the face indices for that cell*/
Eigen::VectorXd Mesh::pyramid_volumes(int cell_nr) const {
  int numFaces = cells[cell_nr].faces.size();
  Eigen::ArrayXi cell_faces = cells[cell_nr].faces;
  Eigen::Vector3d raw_avg = cell_raw_avg(cell_nr);
  Eigen::VectorXd pyr_vols(numFaces);

  Eigen::Vector3d face_cent;
  double pyramid_height;

  for (int i = 0; i < numFaces; i++) {
    face_cent = faces[cell_faces[i]].centroid;
    pyramid_height = (face_cent - raw_avg).norm();

    // Volume of pyramid for which base is polygon is calculated using standard
    // formula
    pyr_vols[i] =
        (1.0 / 3.0) * pyramid_height * (faces[cell_faces[i]].Sf.norm());
  }

  return pyr_vols;
}

// Cell volume is stored as attribute of each cell. Computed
// simply as sum of pyramid volumes for given cell
void Mesh::set_cell_vol(int cell_nr) {
  double cell_volume = 0.0;
  int numFaces = cells[cell_nr].faces.size();
  Eigen::VectorXd pyr_volumes(numFaces);
  pyr_volumes = pyramid_volumes(cell_nr);

  for (int i = 0; i < numFaces; i++) {
    cell_volume += pyr_volumes[i];
  }

  cells[cell_nr].vol = cell_volume;
}

/* Cell centroid will be calculated as the average of face-centroids,
weighted by the volume of the pyramid formed between face and cell raw
average.*/
void Mesh::set_cell_centroid(int cell_nr) {
  int numFaces = cells[cell_nr].faces.size();
  Eigen::VectorXd pyr_volumes = pyramid_volumes(cell_nr);
  Eigen::Vector3d the_centroid(0.0, 0.0, 0.0);
  Eigen::ArrayXi cell_faces;
  cell_faces = cells[cell_nr].faces;
  Eigen::VectorXd f_centroid;
  double cell_volume = 0.0;

  // Compute weighted average
  for (int i = 0; i < numFaces; i++) {
    f_centroid = faces[cell_faces[i]].centroid;
    the_centroid += pyr_volumes[i] * (f_centroid);
    cell_volume += pyr_volumes[i];
  }

  the_centroid /= cell_volume;
  cells[cell_nr].centroid = the_centroid;
}

/*Set face interpolation factor for given face. When interpolating cell variable
ø onto faces, we use formula (f_x)*ø_ownerCell + (1 - f_x)*ø_nbrCell, where f_x
is interpolation factor. Let Nf denote vector between neighbour cell centroid
and face centroid. Let PN denote vector between owner cell centroid and
neighbour cell centroid. We then set f_x = |Nf|/|PN|. We will only invoke this
function for faces not on boundary */
void Mesh::set_face_interpolation_factor(int face_nr) {
  int cellP = faces[face_nr].ownerCell;
  int cellN = faces[face_nr].neighbourCell;

  // Get centroids
  Eigen::Vector3d face_cent = faces[face_nr].centroid;
  Eigen::Vector3d cellP_cent = cells[cellP].centroid;
  Eigen::Vector3d cellN_cent = cells[cellN].centroid;

  // Calculate as outlined above
  double f_x =
      (face_cent - cellN_cent).norm() / (cellP_cent - cellN_cent).norm();
  faces[face_nr].f_x = f_x;
}

/*Function to set the flux across a face*/
void Mesh::calc_face_flux(std::vector<int> b_type,
                          Eigen::MatrixXd boundary_values_v) {
  for (int i = 0; i < faces.size(); i++) {
    if (!faces[i].empty_face) {
      int boundary_id = faces[i].boundary;
      Eigen::Vector3d face_velocity;
      Eigen::Vector3d Sf = faces[i].Sf;
      int owner = faces[i].ownerCell;

      // If face is on boundary, determine boundary condition
      // type and set flux accordingly
      if (boundary_id >= 0) {
        if (b_type[boundary_id] == 0) {
          // Dirichlet Boundary
          face_velocity = boundary_values_v.col(boundary_id);
        }

        else {
          // Neumann Boundary
          face_velocity =
              (1 / faces[i].delta) * boundary_values_v.col(boundary_id) +
              cells[owner].vel;
        }

        faces[i].flux = face_velocity.dot(Sf);
      }

      /* else face is internal. Interpolate velocities
      from cell centres to face*/
      else {
        int nbr = faces[i].neighbourCell;
        const double f_x = faces[i].f_x;
        face_velocity = f_x * cells[owner].vel + (1 - f_x) * cells[nbr].vel;

        faces[i].flux = Sf.dot(face_velocity);
      }
    }
  }
}

/*Function to correct face flux*/
void Mesh::correct_face_flux(Eigen::VectorXd scaler, std::vector<int> b_type,
                             std::vector<double> boundary_values_p) {
  for (int i = 0; i < faces.size(); i++) {
    if (!faces[i].empty_face) {
      Face& face = faces[i];
      int boundary_id = face.boundary;
      double sf_grad_p;
      double scaler_face;
      Eigen::Vector3d Sf = face.Sf;
      int owner = face.ownerCell;

      // If face is on boundary, determine boundary condition
      // type and correct fluc accordingly using pressure b-value
      if (boundary_id >= 0) {
        if (b_type[boundary_id] == 0) {
          // Dirichlet Boundary on pressure
          // coef is estimate of dot product of face normal with pressure
          // gradient
          sf_grad_p = Sf.norm() *
                      (boundary_values_p[boundary_id] - cells[owner].p) *
                      face.delta;

          face.flux -= scaler(owner) * sf_grad_p;
        }

        else {
          // Neumann Boundary
          sf_grad_p = Sf.norm() * boundary_values_p[boundary_id];
          face.flux -= scaler(owner) * sf_grad_p;
        }
      }

      /* else face is internal. Interpolate pressure gradient
      from pressure at adjacent cell centres and account for potential
      non-orthogonalit*/
      else {
        int nbr = face.neighbourCell;
        double f_x = face.f_x;
        double coef = Sf.norm() * (cells[nbr].p - cells[owner].p) * face.delta;
        double scaler_face = f_x * scaler(owner) + (1 - f_x) * scaler(nbr);

        face.flux -= scaler_face * coef;

        // Account for potential orthogonal correction
        //  Compute vector between cell centroids
        // Eigen::Vector3d PN = cells[nbr].centroid - cells[owner].centroid;

        // Compute correction k_f
        // double alpha = Sf.dot(Sf) / PN.dot(Sf);
        // Eigen::Vector3d k = Sf - alpha * PN;

        // Due to numerical accuracy, set lower bound for |k_f|,
        // below which we force k_f to zero
        // if (k.norm() < 1e-10) {
        //  k = Eigen::Vector3d::Zero(3);
        //}

        // Account for non-orthogonal correction
        // Only old gradients are involved in this calculation
        // so we subtract from b-vector (source term)
        // Eigen::Vector3d gradp_face_old =
        //   f_x * cells[owner].grad_p + (1 - f_x) * cells[nbr].grad_p;

        // sf_grad_p =
        //(alpha * PN.norm()) * face.delta * (cells[nbr].p - cells[owner].p) +
        // k.dot(gradp_face_old);
      }
    }
  }
}

double Mesh::get_face_flux(int face_nr) const { return faces[face_nr].flux; }

Eigen::Vector3d Mesh::get_face_centroid(int face_nr) const {
  return faces[face_nr].centroid;
}

void Mesh::correct_velocity(Eigen::VectorXd scaler) {
  for (int i = 0; i < cells.size(); i++) {
    cells[i].vel -= scaler(i) * cells[i].grad_p;
  }
}

/*Set delta coefficient for given face. If face is not on boundary,
delta coefficient is simply 1/(distance between adjoining cell centres).
Otherwise, delta coefficient is 1/(distance from owner cell centroid to face
centroid)*/
void Mesh::set_delta_coef(int face_nr) {
  int cellP = faces[face_nr].ownerCell;
  int cellN = faces[face_nr].neighbourCell;

  // If face is on boundary (has no neighbour cell)
  if (cellN < 0) {
    Eigen::Vector3d cellP_cent = cells[cellP].centroid;
    Eigen::Vector3d face_cent = faces[face_nr].centroid;
    double delta = 1 / (cellP_cent - face_cent).norm();
    faces[face_nr].delta = delta;
  }

  // If face is not on boundary of domain. Note, empty faces
  // included here, although their delta coefficient will not be used
  else {
    Eigen::Vector3d cellP_cent = cells[cellP].centroid;
    Eigen::Vector3d cellN_cent = cells[cellN].centroid;
    double delta = 1 / (cellP_cent - cellN_cent).norm();
    faces[face_nr].delta = delta;
  }
}

// Calculate the integral over each cell volume of the gradient of a
//  a field. Pass arguments with boundary conditions for the given field
//  at each boundary patch along with boundary condition type
void Mesh::calc_gauss_gradients(std::vector<int>& bc_type,
                                std::vector<double>& boundary_values,
                                Field field) {
  int nr_cells = getNumCells();
  // Loop over cells to calculate gradient at each centroid
  for (int i = 0; i < nr_cells; i++) {
    double val_i = get_cell_val(i, field);

    Eigen::Vector3d gradient = Eigen::Vector3d::Zero(3);

    // Extract the list of faces for a cell
    Eigen::ArrayXi cell_faces = cells[i].faces;

    // Loop over faces
    for (int j = 0; j < cell_faces.size(); j++) {
      if (!faces[cell_faces[j]].empty_face) {
        const Face& face = faces[cell_faces[j]];

        Eigen::Vector3d Sf = face.Sf;

        int boundary = face_boundary(cell_faces[j]);

        // If face not a boundary, interpolate field value to face
        if (boundary < 0) {
          int N = (i == face.ownerCell) ? face.neighbourCell : face.ownerCell;
          double f_x = (i == face.ownerCell) ? face.f_x : (1 - face.f_x);
          double val_N = get_cell_val(N, field);

          if (i == face.neighbourCell) {
            Sf *= -1;
          }

          gradient += (f_x * val_i + (1 - f_x) * val_N) * Sf;

        }

        else {
          // If face at boundary with Dirichlet condition
          if (bc_type[boundary] == 0) {
            gradient += boundary_values[boundary] * Sf;
          }

          // if face at boundary with Neumann condition
          else {
            double dist = 1 / face.delta;
            double val_f = val_i + dist * boundary_values[boundary];
            gradient += val_f * Sf;
          }
        }
      }
    }
    set_cell_grad(i, field, gradient);
  }
}
