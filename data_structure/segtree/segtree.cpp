template <typename T, auto op, auto e>
struct SegTree {
   public:
    SegTree(int n) : n(n), lgn(log2_ceil(n)), sz(1 << lgn), tree(sz << 1, e()) {
        init();
    }
    SegTree(vector<T>& v)
        : n(v.size()), lgn(log2_ceil(n)), sz(1 << lgn), tree(sz << 1, e()) {
        for (int i = 0; i < (int)v.size(); i++) tree[i | sz] = v[i];
        init();
    }
    void update(int i, const T& val) {
        tree[i |= sz] = val;
        while (i >>= 1) tree[i] = op(tree[i << 1], tree[i << 1 | 1]);
    }
    T point_query(int idx) const { return tree[idx | sz]; }
    T range_query(int l, int r) const {
        T res = e();
        for (l |= sz, r |= sz; l <= r; l >>= 1, r >>= 1) {
            if (l & 1) res = op(res, tree[l++]);
            if (~r & 1) res = op(res, tree[r--]);
        }
        return res;
    }

   private:
    const int n, lgn, sz;
    vector<T> tree;
    static int log2_ceil(int n) {
        int res = 0;
        while (n > 1 << res) res++;
        return res;
    }
    void init() {
        for (int i = sz - 1; i; i--)
            tree[i] = op(tree[i << 1], tree[i << 1 | 1]);
    }
};