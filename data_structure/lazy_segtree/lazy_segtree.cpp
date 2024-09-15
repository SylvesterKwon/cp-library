template <typename T, auto op, auto e, typename F, auto apply_F,
          auto F_composition, auto F_id>
struct LazySegTree {
   public:
    LazySegTree(vector<T>& v)
        : n(v.size()), lgn(log2_ceil(n)), sz(1 << lgn), tree(sz << 1, e()) {
        for (int i = 0; i < n; i++) tree[i | sz] = v[i];
        init();
    }
    void set(int idx, const T& value) {  // need qa
        idx |= sz;
        for (int i = lgn; i; i--) flush_lazy(idx >> i);
        tree[idx] = value;
        for (int i = 1; i <= lgn; i++) apply_op(idx >> i);
    }
    void range_update(int l, int r, const F& f) {
        l |= sz, r |= sz;
        for (int i = lgn; i; i--) {
            if (l >> i << i != l) flush_lazy(l >> i);
            if ((r + 1) >> i << i != r + 1) flush_lazy(r >> i);
        }
        for (int ll = l, rr = r; ll <= rr; ll >>= 1, rr >>= 1) {
            if (ll & 1) evalutate_update(ll++, f);
            if (~rr & 1) evalutate_update(rr--, f);
        }
        for (int i = 1; i <= lgn; i++) {
            if (l >> i << i != l) apply_op(l >> i);
            if ((r + 1) >> i << i != r + 1) apply_op(r >> i);
        }
    }
    T point_query(int idx) {  // need qa
        idx |= sz;
        for (int i = lgn; i; i--) flush_lazy(idx >> i);
        return tree[idx];
    }
    T range_query(int l, int r) {
        T res = e();
        l |= sz, r |= sz;
        for (int i = lgn; i; i--) {
            if (l >> i << i != l) flush_lazy(l >> i);
            if ((r + 1) >> i << i != r + 1) flush_lazy(r >> i);
        }
        for (int ll = l, rr = r; ll <= rr; ll >>= 1, rr >>= 1) {
            if (ll & 1) res = op(res, tree[ll++]);
            if (~rr & 1) res = op(res, tree[rr--]);
        }
        return res;
    }

   private:
    const int n, lgn, sz;
    vector<T> tree;
    vector<F> lazy_tree;
    static int log2_ceil(int n) {
        int res = 0;
        while (n > 1 << res) res++;
        return res;
    }
    void init() {
        lazy_tree.resize(sz, F_id());
        for (int i = sz - 1; i; i--) apply_op(i);
    }
    void evalutate_update(int idx, const F& f) {
        tree[idx] = apply_F(f, tree[idx]);
        if (idx < sz) lazy_tree[idx] = F_composition(f, lazy_tree[idx]);
    }
    void flush_lazy(int idx) {
        evalutate_update(idx << 1, lazy_tree[idx]);
        evalutate_update(idx << 1 | 1, lazy_tree[idx]);
        lazy_tree[idx] = F_id();
    }
    void apply_op(int idx) {
        tree[idx] = op(tree[idx << 1], tree[idx << 1 | 1]);
    }
};