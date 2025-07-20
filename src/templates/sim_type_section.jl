using Mustache

types_tmp = mt"""
#ifndef DGGML_PLANT_TYPES_HPP
#define DGGML_PLANT_TYPES_HPP

#include "YAGL_Graph.hpp"
#include "YAGL_Node.hpp"
#include "SpatialData3D.hpp"

namespace Plant
{
    const uint8_t DIM3D = 3;

    struct Negative
    {
        double velocity[DIM3D];
        double unit_vec[DIM3D];
    };

    struct Intermediate
    {
        double velocity[DIM3D];
        double unit_vec[DIM3D];
    };

    struct Positive
    {
        double velocity[DIM3D];
        double unit_vec[DIM3D];
    };

    struct Junction
    {
        double unit_vec[DIM3D];
    };

    struct Holding
    {
        double unit_vec[DIM3D];
        double collision_angle;
    };

    struct Capture
    {
        double unit_vec[DIM3D];
    };

    struct Zipper
    {
        double unit_vec[DIM3D];
    };

    struct Boundary {};

    struct Nucleator {};

    using key_type = std::size_t;

    using graph_type = YAGL::Graph<key_type,
    SpatialNode3D<Negative, Intermediate, Positive, Junction, Zipper, Holding, Capture, Boundary, Nucleator>>;
}

#endif
"""

function generate_type()
    items = Dict(
             "vec" => ["a\n", "b\n", "c\n"]
            )
    
    # res = Mustache.render(section, items)
    # println("$res")
    # res
    tpl = mt"oof {{#:vec}}{{.}}{{{\n}}}{{/:vec}} oof" # just first one
    res = Mustache.render(tpl, items)
    println("$res")
end

function generate_types_section(data)
    generate_type()
    Mustache.render(types_tmp, data)
end

generate_type()
