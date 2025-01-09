// Copyright 2011-2020 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#include "pmp/algorithms/sec.h"

#include <iterator>
#include <limits>
#include <numbers>

#include "pmp/algorithms/distance_point_triangle.h"
#include "pmp/algorithms/normals.h"

namespace pmp {
namespace {

template <class HeapEntry, class HeapInterface>
class Heap : private std::vector<HeapEntry>
{
public:
    using This = Heap<HeapEntry, HeapInterface>;

    // Constructor
    Heap() : HeapVector() {}

    // Construct with a given \p HeapInterface.
    Heap(const HeapInterface& i) : HeapVector(), interface_(i) {}

    // Destructor.
    ~Heap() = default;

    // clear the heap
    void clear() { HeapVector::clear(); }

    // is heap empty?
    bool empty() { return HeapVector::empty(); }

    // returns the size of heap
    unsigned int size() { return (unsigned int)HeapVector::size(); }

    // reserve space for N entries
    void reserve(unsigned int n) { HeapVector::reserve(n); }

    // reset heap position to -1 (not in heap)
    void reset_heap_position(HeapEntry h)
    {
        interface_.set_heap_position(h, -1);
    }

    // is an entry in the heap?
    bool is_stored(HeapEntry h)
    {
        return interface_.get_heap_position(h) != -1;
    }

    // insert the entry h
    void insert(HeapEntry h)
    {
        This::push_back(h);
        upheap(size() - 1);
    }

    // get the first entry
    HeapEntry front()
    {
        assert(!empty());
        return entry(0);
    }

    // delete the first entry
    void pop_front()
    {
        assert(!empty());
        interface_.set_heap_position(entry(0), -1);
        if (size() > 1)
        {
            entry(0, entry(size() - 1));
            HeapVector::resize(size() - 1);
            downheap(0);
        }
        else
            HeapVector::resize(size() - 1);
    }

    // remove an entry
    void remove(HeapEntry h)
    {
        int pos = interface_.get_heap_position(h);
        interface_.set_heap_position(h, -1);
        if (pos == -1)
        {
            return;
        }
        assert(pos != -1);
        assert((unsigned int)pos < size());

        // last item ?
        if ((unsigned int)pos == size() - 1)
            HeapVector::resize(size() - 1);

        else
        {
            entry(pos, entry(size() - 1)); // move last elem to pos
            HeapVector::resize(size() - 1);
            downheap(pos);
            upheap(pos);
        }
    }

    // update an entry: change the key and update the position to
    // reestablish the heap property.
    void update(HeapEntry h)
    {
        int pos = interface_.get_heap_position(h);
        assert(pos != -1);
        assert((unsigned int)pos < size());
        downheap(pos);
        upheap(pos);
    }

    // Check heap condition. true if heap condition is satisfied, false if not.
    bool check()
    {
        bool ok(true);
        unsigned int i, j;
        for (i = 0; i < size(); ++i)
        {
            if (((j = left(i)) < size()) &&
                interface_.greater(entry(i), entry(j)))
            {
                ok = false;
            }
            if (((j = right(i)) < size()) &&
                interface_.greater(entry(i), entry(j)))
            {
                ok = false;
            }
        }
        return ok;
    }

private:
    using HeapVector = std::vector<HeapEntry>;

    // Upheap. Establish heap property.
    void upheap(unsigned int idx)
    {
        HeapEntry h = entry(idx);
        unsigned int parentIdx;

        while ((idx > 0) && interface_.less(h, entry(parentIdx = parent(idx))))
        {
            entry(idx, entry(parentIdx));
            idx = parentIdx;
        }

        entry(idx, h);
    }

    // Downheap. Establish heap property.
    void downheap(unsigned int idx)
    {
        HeapEntry h = entry(idx);
        unsigned int childIdx;
        unsigned int s = size();

        while (idx < s)
        {
            childIdx = left(idx);
            if (childIdx >= s)
                break;

            if ((childIdx + 1 < s) &&
                (interface_.less(entry(childIdx + 1), entry(childIdx))))
                ++childIdx;

            if (interface_.less(h, entry(childIdx)))
                break;

            entry(idx, entry(childIdx));
            idx = childIdx;
        }

        entry(idx, h);
    }

    // Get the entry at index idx
    inline HeapEntry entry(unsigned int idx)
    {
        assert(idx < size());
        return (This::operator[](idx));
    }

    // Set entry H to index idx and update H's heap position.
    inline void entry(unsigned int idx, HeapEntry h)
    {
        assert(idx < size());
        This::operator[](idx) = h;
        interface_.set_heap_position(h, idx);
    }

    // Get parent's index
    inline unsigned int parent(unsigned int i) { return (i - 1) >> 1; }

    // Get left child's index
    inline unsigned int left(unsigned int i) { return (i << 1) + 1; }

    // Get right child's index
    inline unsigned int right(unsigned int i) { return (i << 1) + 2; }

    // Instance of HeapInterface
    HeapInterface interface_;
};

class sec
{
public:
    sec(SurfaceMesh& mesh);
    ~sec();

    void decimate(unsigned int n_vertices);

private:
    // Store data for an halfedge collapse
    struct CollapseData
    {
        CollapseData(SurfaceMesh& sm, Halfedge h);

        SurfaceMesh& mesh;

        /*        vl
         *        *
         *       / \
         *      /   \
         *     / fl  \
         * v0 *------>* v1
         *     \ fr  /
         *      \   /
         *       \ /
         *        *
         *        vr
         */
        Halfedge h;
        Halfedge v0v1; // Halfedge to be collapsed
        Halfedge v1v0; // Reverse halfedge
        Vertex v0;     // Vertex to be removed
        Vertex v1;     // Remaining vertex
        Face fl;       // Left face
        Face fr;       // Right face
        Vertex vl;     // Left vertex
        Vertex vr;     // Right vertex
        Halfedge v1vl, vlv0, v0vr, vrv1;
    };

    // Heap interface
    class HeapInterface
    {
    public:
        HeapInterface(EdgeProperty<float> prio, EdgeProperty<int> pos)
            : prio_(prio), pos_(pos)
        {
        }

        bool less(Edge e0, Edge e1) { return prio_[e0] < prio_[e1]; }
        bool greater(Edge e0, Edge e1) { return prio_[e0] > prio_[e1]; }
        int get_heap_position(Edge e) { return pos_[e]; }
        void set_heap_position(Edge e, int pos) { pos_[e] = pos; }

    private:
        EdgeProperty<float> prio_;
        EdgeProperty<int> pos_;
    };

    using PriorityQueue = Heap<Edge, HeapInterface>;

    using Points = std::vector<Point>;

    // put the vertex v in the priority queue
    void enqueue_edge(PriorityQueue& queue, Edge e);

    // is collapsing the halfedge h allowed?
    bool is_collapse_legal(const CollapseData& cd);

    float priority(const Edge e);
    float priority(const CollapseData& cd);

    SurfaceMesh& mesh_;

    EdgeProperty<float> epriority_;
    EdgeProperty<int> heap_pos_;

    VertexProperty<Point> vpoint_;
};

sec::sec(SurfaceMesh& mesh) : mesh_(mesh)
{
    if (!mesh_.is_triangle_mesh())
        throw InvalidInputException("Input is not a triangle mesh!");

    // get properties
    vpoint_ = mesh_.vertex_property<Point>("v:point");
}

sec::~sec() {}

void sec::decimate(unsigned int n_vertices)
{
    std::vector<Edge> one_ring;

    epriority_ = mesh_.add_edge_property<float>("e:prio");
    heap_pos_ = mesh_.add_edge_property<int>("e:heap");

    HeapInterface hi(epriority_, heap_pos_);
    PriorityQueue queue(hi);
    queue.reserve(mesh_.n_edges());
    for (auto e : mesh_.edges())
    {
        queue.reset_heap_position(e);
        enqueue_edge(queue, e);
    }

    auto nv = mesh_.n_vertices();
    while (nv > n_vertices && !queue.empty())
    {
        auto e = queue.front();
        queue.pop_front();

        if (mesh_.is_deleted(e))
        {
            continue;
        }

        Halfedge h = mesh_.halfedge(e, 0);

        assert(mesh_.edge(h) == e);

        assert(mesh_.is_valid(h));
        assert(!mesh_.is_deleted(h));

        CollapseData cd(mesh_, h);

        if (!mesh_.is_collapse_ok(h))
        {
            continue;
        }

        auto mid_point = (mesh_.vertex_property<Point>("v:point")[cd.v0] +
                          mesh_.vertex_property<Point>("v:point")[cd.v1]) /
                         2.0;

        one_ring.clear();
        for (auto ee : mesh_.edges(cd.v0))
        {
            if (ee == e)
            {
                continue;
            }
            one_ring.push_back(ee);
        }

        for (auto ee : mesh_.edges(cd.v1))
        {
            if (ee == e)
            {
                continue;
            }
            one_ring.push_back(ee);
        }

        mesh_.collapse(h);
        --nv;

        assert(!mesh_.is_deleted(cd.v1));

        for (auto ee : mesh_.edges(cd.v1))
        {
            if (ee == e)
            {
                continue;
            }
            one_ring.push_back(ee);
        }

        mesh_.position(cd.v1) = mid_point;

        for (auto ee : one_ring)
        {
            if (!mesh_.is_deleted(ee))
            {
                enqueue_edge(queue, ee);
            }
            else
            {
                queue.remove(ee);
            }
        }
    }

    mesh_.garbage_collection();
    mesh_.remove_edge_property(epriority_);
    mesh_.remove_edge_property(heap_pos_);
}

void sec::enqueue_edge(PriorityQueue& queue, Edge e)
{
    std::vector<Halfedge> he;

    assert(mesh_.is_valid(e));
    assert(!mesh_.is_deleted(e));

    Halfedge h0 = mesh_.halfedge(e, 0);
    Halfedge h1 = mesh_.halfedge(e, 1);

    he.push_back(h0);
    he.push_back(h1);

    bool is_legal = true;

    for (auto h : he)
    {
        assert(mesh_.is_valid(h));
        assert(!mesh_.is_deleted(h));

        CollapseData cd(mesh_, h);
        if (!is_collapse_legal(cd))
        {
            is_legal = false;
        }
    }

    if (is_legal)
    {
        epriority_[e] = priority(e);

        if (queue.is_stored(e))
        {
            queue.update(e);
        }
        else
        {
            queue.insert(e);
        }
    }
}

bool sec::is_collapse_legal(const CollapseData& cd)
{
    // do not collapse boundary vertices to interior vertices
    if (mesh_.is_boundary(cd.v0) && !mesh_.is_boundary(cd.v1))
        return false;

    // there should be at least 2 incident faces at v0
    if (mesh_.cw_rotated_halfedge(mesh_.cw_rotated_halfedge(cd.v0v1)) ==
        cd.v0v1)
        return false;

    // topological check
    if (!mesh_.is_collapse_ok(cd.v0v1))
        return false;

    // collapse passed all tests -> ok
    return true;
}

float sec::priority(const CollapseData& cd)
{
    Edge e = cd.mesh.edge(cd.h);

    return priority(e);
}

float sec::priority(const Edge e)
{
    assert(mesh_.is_valid(e));
    assert(!mesh_.is_deleted(e));

    auto v0 = mesh_.vertex(e, 0);
    auto v1 = mesh_.vertex(e, 1);

    assert(mesh_.is_valid(v0));
    assert(!mesh_.is_deleted(v0));

    assert(mesh_.is_valid(v1));
    assert(!mesh_.is_deleted(v1));

    float l = distance(mesh_.vertex_property<Point>("v:point")[v0],
                       mesh_.vertex_property<Point>("v:point")[v1]);
    return l;
}

sec::CollapseData::CollapseData(SurfaceMesh& sm, Halfedge h) : mesh(sm), h(h)
{
    v0v1 = h;
    v1v0 = mesh.opposite_halfedge(v0v1);
    v0 = mesh.to_vertex(v1v0);
    v1 = mesh.to_vertex(v0v1);
    fl = mesh.face(v0v1);
    fr = mesh.face(v1v0);

    // get vl
    if (fl.is_valid())
    {
        v1vl = mesh.next_halfedge(v0v1);
        vlv0 = mesh.next_halfedge(v1vl);
        vl = mesh.to_vertex(v1vl);
    }

    // get vr
    if (fr.is_valid())
    {
        v0vr = mesh.next_halfedge(v1v0);
        vrv1 = mesh.prev_halfedge(v0vr);
        vr = mesh.from_vertex(vrv1);
    }
}
} // namespace

void shortest_edge_collapse(SurfaceMesh& mesh, unsigned int n_vertices)
{
    sec decimator(mesh);

    decimator.decimate(n_vertices);
}

} // namespace pmp
