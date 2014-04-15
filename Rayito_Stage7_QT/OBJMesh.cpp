#include <fstream>
#include <sstream>
#include <iostream>

#include <RMesh.h>


namespace Rayito
{

/*
 * Brief overview of OBJ file format.  It is an ASCII format.  It has comments:
 *     # This is a comment
 * Vertices:
 *     v 0.1 -2.0 3.5
 *   -or-
 *     v 0.1 -2.0 3.5 1.0  # Homogeneous coordinates, rarely encountered in OBJ
 * Normals:
 *     vn 1.0 0.0 0.0
 * UVs / texture coordinates:
 *     vt 0.1 0.2
 *   -or-
 *     vt 0.1 0.2 0.3  # 3D or homogeneous UVs, rarely encountered in OBJ
 * Faces:
 *     f 0 1 2 3
 *   -or-
 *     f 0/0/0 1/2/3 4/4/5
 *   -or-
 *     f 0//0 1//3 4//5
 *   -or-
 *     f 0/0 1/2 4/4
 * Form 1 is just vertex indices, form 2 is vertex/UV/normal indices, form 3 is
 * vertex and normal indices, and form 4 is vertex and UV indices.  The indices
 * are 1-based, so to tranlate to C arrays we have to subtract 1 from each index.
 * If any index is missing (normal and/or UV) then those have to be missing from
 * ALL faces, they have to all be the same style of spec.  Indices can also be
 * negative, in which case they index backwards from the last vertex/UV/normal
 * that was listed.
 * 
 * There are various other directives which we ignore; they are used for listing
 * multiple objects in the file (o, s), and materials (mtllib, usemtl) which
 * refer to a .mtl companion file.  Implementing those are left as an exercise
 * to the reader.
 * 
 * Also note that for this stage, we do not support texture mapping, so the vt
 * directive is also effectively ignored.  But the code is in there, so it
 * should be easy to enable it in the future.
 */
Mesh* createFromOBJFile(const char* filename)
{
    std::ifstream input(filename);
    
    std::vector<Face> faces;
    std::vector<Point> verts;
    std::vector<Vector> normals;
    
    std::string lineStr;
    std::string command;
    while (input.good())
    {
        lineStr.clear();
        std::getline(input, lineStr);

        std::istringstream lineInput(lineStr);
        if (lineInput.eof())
        {
            continue;
        }

        command.clear();
        lineInput >> command;
        if (lineInput.fail())
        {
            continue;
        }

        if (command[0] == '#')
        {
            // Found a comment; eat it
        }
        else if (command == "v")
        {
            // NOTE: there is an optional w coordinate that we're ignoring here
            Point v;
            lineInput >> v.m_x;
            lineInput >> v.m_y;
            lineInput >> v.m_z;
            verts.push_back(v);
        }
        else if (command == "vn")
        {
            Vector v;
            lineInput >> v.m_x;
            lineInput >> v.m_y;
            lineInput >> v.m_z;
            normals.push_back(v);
        }
        else if (command == "vt")
        {
            // Note: there's an optional w coordinate that we're ignoring here
//            Vec2 uv;
//            lineInput >> uv.m_u;
//            lineInput >> uv.m_v;
//            uvs.push_back(uv);
        }
        else if (command == "f")
        {
            faces.push_back(Face());
            Face& face = faces.back();
            while (lineInput.good())
            {
                int vi;
                lineInput >> vi;
                if (lineInput.fail())
                    break;
                int uvi, ni;
                bool gotUV = false;
                bool gotN = false;
                if (lineInput.peek() == '/')
                {
                    char slash;
                    lineInput >> slash;
                    if (lineInput.peek() == '/')
                    {
                        lineInput >> slash;
                        lineInput >> ni;
                        gotN = true;
                    }
                    else
                    {
                        lineInput >> uvi;
                        gotUV = true;
                        if (lineInput.peek() == '/')
                        {
                            lineInput >> slash;
                            lineInput >> ni;
                            gotN = true;
                        }
                    }
                }
                vi = vi > 0 ? vi - 1 : (int)verts.size() + vi;
                face.m_vertexIndices.push_back(vi);
                if (vi >= (int)verts.size())
                    std::cerr << "Found out-of-range vertex index: " << vi << std::endl;
                if (gotUV)
                {
//                    uvi = uvi > 0 ? uvi - 1 : (int)uvs.size() + uvi;
//                    face.m_uvIndices.push_back(uvi);
//                    if (uvi >= uvs.size())
//                        std::cerr << "Found out-of-range UV index: " << uvi << std::endl;
                }
                if (gotN)
                {
                    ni = ni > 0 ? ni - 1 : (int)normals.size() + ni;
                    face.m_normalIndices.push_back(ni);
                    if (ni >= (int)normals.size())
                        std::cerr << "Found out-of-range N index: " << ni << std::endl;
                }
            }
        }
        else if (command == "usemtl")
        {
            
        }
        else if (command == "mtllib")
        {
            
        }
        else if (command == "s")
        {
            
        }
        else
        {
            
        }
    }
    if (verts.empty() || faces.empty())
        return NULL;
    return new Mesh(verts, normals, faces, NULL);
}


} // namespace Rayito
