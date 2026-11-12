"""Extraction on the e-graph read back from egglog. Each e-class's minimal-cost set F_k(c) is a
list of (cost, program) entries whose costs are Min(X_k(c)).

A round keeps Min(F_k(c) | X_{k+1}(c)), which is Min(X_{k+1}(c)) for a monotone cost. A cost that is
FAIL is left out of every minimal-cost set."""

import time
from itertools import product

from .analysis import FAIL, err, precedes


def insert(front, entry, leq):
    """Add (cost, program) to a set of minimal costs unless a cost there is <= its cost, and drop
    the costs it is below. A cost keeps the first program found for it."""
    if not any(leq(e, entry[0]) for e, _ in front):
        front[:] = [f for f in front if not leq(entry[0], f[0])] + [entry]


def program(node, programs):
    """The program of a node, from its children's."""
    match node.op:
        case "Var":
            return ("Var", node.payload)
        case "Num" | "Lit":
            return ("Num", node.payload)
        case "Ctx":
            return programs[0]
        case "If":
            return ("If", (">", *node.payload), *programs)
        case op:
            return (op, *programs)


def assumes(node, box):
    """The subdomain of an e-node, as the set of guard sides (x, k, side) it assumes, where side
    is True for x > k and False for x <= k."""
    if any(child not in box for child in node.children):
        return None
    guards = [box[child] for child in node.children]
    match node.op:
        case "Ctx":
            return guards[0] | {node.payload}
        case "If":
            x, k = node.payload
            return (guards[0] - {(x, k, True)}) | (guards[1] - {(x, k, False)})
        case _:
            return frozenset().union(*guards)


def boxes(g):
    """Each e-class's subdomain V_c, as a set of guard sides: the intersection of its e-nodes'."""
    box = {}
    while True:
        new = {}
        for c, nodes in g.nodes.items():
            guards = [s for s in (assumes(n, box) for n in nodes) if s is not None]
            if guards:
                new[c] = frozenset.intersection(*guards)
        if new == box:
            return box
        box = new


def participants(g):
    """The e-nodes whose subdomain is their e-class's V_c, so that V_c lies inside each child's."""
    box = boxes(g)
    nodes = [(c, n) for c, ns in g.nodes.items() for n in ns
             if (guards := assumes(n, box)) is not None and guards <= box[c]]
    members, readers = {c: [] for c in g.nodes}, {c: set() for c in g.nodes}
    for i, (c, n) in enumerate(nodes):
        members[c].append(i)
        for child in n.children:
            readers[child].add(i)
    return nodes, members, readers


def extract(g, time_limit=None, cost=err, leq=precedes):
    """F_k of every e-class for the largest k reached, k, and whether the time limit stopped it
    before the fixed point F_{k+1} = F_k. Round k+1 evaluates again only the e-nodes with a child
    whose minimal-cost set changed in round k."""
    deadline = time_limit and time.monotonic() + time_limit
    nodes, members, readers = participants(g)

    def costs(c, node):
        """Min of what node adds to X_{k+1}(c), from its children's F_k, or None at the deadline."""
        Xs = [g.interval[child] for child in node.children]
        out = []
        for i, entries in enumerate(product(*(F[child] for child in node.children))):
            if i % 1024 == 0 and deadline and time.monotonic() > deadline:
                return None
            e = cost(node.op, g.interval[c], Xs, [e for e, _ in entries])
            if e is not FAIL:
                insert(out, (e, program(node, [p for _, p in entries])), leq)
        return out

    F = {c: [] for c in g.nodes}  # F_0
    added, todo, k = {}, range(len(nodes)), 0
    while todo:
        for i in todo:
            added[i] = costs(*nodes[i])
            if added[i] is None:
                return F, k, True
        changed = {}
        for c in {nodes[i][0] for i in todo}:
            front = list(F[c])
            for i in members[c]:
                for entry in added[i]:
                    insert(front, entry, leq)
            if {e for e, _ in front} != {e for e, _ in F[c]}:
                changed[c] = front
        F.update(changed)
        k += 1
        todo = sorted({i for c in changed for i in readers[c]})
    return F, k, False


def best(front, bound):
    """The entry minimizing a monotone readout, such as a bound from analysis.BOUNDS."""
    return min(front, key=lambda entry: bound(entry[0]))


def analyze(g, p):
    """The cost of one program p as a term of g, with its ifs written as egraph.guarded writes them.
    Each e-class's range is the one the e-graph gives it, so the bound can be tighter than the
    program's own analysis: the ranges of the other terms in its e-classes narrow it."""
    index = {(n.op, n.payload, n.children): c for c, ns in g.nodes.items() for n in ns}

    def node(op, payload, kids):
        c = index[(op, payload, tuple(k for k, _ in kids))]
        Xs = [g.interval[k] for k, _ in kids]
        return c, err(op, g.interval[c], Xs, [e for _, e in kids])

    def walk(p, guard):
        match p:
            case ("Var", x) if x in guard:
                return node("Ctx", guard[x], [walk(p, {})])
            case ("Var" | "Num", v):
                return node(p[0] if (p[0], v, ()) in index else "Lit", v, [])
            case ("If", (_, x, k), a, b):
                return node("If", (x, k), [walk(a, guard | {x: (x, k, True)}), walk(b, guard | {x: (x, k, False)})])
            case (op, *args):
                return node(op, None, [walk(a, guard) for a in args])

    return walk(p, {})[1]
