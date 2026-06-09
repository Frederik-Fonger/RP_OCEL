from collections import Counter
import networkx as nx
from tqdm import tqdm
from multiprocessing import Pool
import networkx as nx
from tqdm import tqdm
import math

def calculate_preset(eog):
    preset = {}
    for e in tqdm(eog.nodes):
        #several different ways of calculating
        # USE THIS FOR LARGE EVENT LOGS
        preset[e] = list(nx.ancestors(eog,e))


        #stable speed also for later events, large logs with large connected components
        #preset[e] = [v for v in nx.dfs_predecessors(EOG, source=e).keys() if v!=e]

        #fast for small graphs/no connected component
        #preset[e] = [n for n in nx.traversal.bfs_tree(EOG, e, reverse=True) if n != e]
    return preset


def calculate_ancestors(args):
    node, graph = args
    return node, list(nx.ancestors(graph, node))


def calculate_preset_mp(eog):
    preset = {}

    # Erstelle die Argumentliste für den Pool

    args = [(node, eog) for node in eog.nodes]
    # Verwende einen Multiprocessing-Pool
    with Pool(processes=5) as pool:
        # Verwende imap für einen Progress-Bar mit tqdm
        results = list(tqdm(pool.imap(calculate_ancestors, args), total=len(args)))

    # Konvertiere die Ergebnisse in ein Dictionary
    preset = dict(results)

    return preset


def calculate_contexts_and_bindings(ocel):
    log = ocel.log.log.copy()
    object_types = ocel.object_types
    contexts = {}
    bindings = {}
    preset = calculate_preset(ocel.graph.eog)
    log["event_objects"] = log.apply(lambda x: [(ot,o) for ot in object_types for o in x[ot]], axis = 1)
    exploded_log = log.explode("event_objects")
    counter_e=0
    for event in tqdm(preset.keys()):
        counter_e+=1
        context = {}
        obs = list(set().union(*log.loc[log["event_id"].isin(preset[event]+[event])]["event_objects"].to_list()))
        binding_sequence = log.loc[log["event_id"].isin(preset[event])].apply(lambda y: (y["event_activity"], { ot : [o for (ot_,o) in y["event_objects"] if ot_ == ot] for ot in object_types}), axis = 1).values.tolist()
        for ob in obs:
            prefix = tuple(exploded_log[(exploded_log["event_objects"] == ob) & (exploded_log["event_id"].isin(preset[event]))]["event_activity"].to_list())
            if ob[0] not in context.keys():
                context[ob[0]] = Counter()
            context[ob[0]]+=Counter([prefix])
        contexts[event] = context
        bindings[event] = binding_sequence
    return contexts, bindings


def process_event(args):
    event, preset_event, log, exploded_log, object_types = args
    context = {}
    obs = list(set().union(*log.loc[log["event_id"].isin(preset_event + [event])]["event_objects"].to_list()))
    binding_sequence = log.loc[log["event_id"].isin(preset_event)].apply(
        lambda y: (y["event_activity"],
                   {ot: [o for (ot_, o) in y["event_objects"] if ot_ == ot] for ot in object_types}
                   ), axis=1).values.tolist()

    for ob in obs:
        prefix = tuple(exploded_log[
                           (exploded_log["event_objects"] == ob) &
                           (exploded_log["event_id"].isin(preset_event))
                           ]["event_activity"].to_list())
        if ob[0] not in context:
            context[ob[0]] = Counter()
        context[ob[0]] += Counter([prefix])

    return event, context, binding_sequence


def calculate_contexts_and_bindings_mp(ocel):
    log = ocel.log.log.copy()
    object_types = ocel.object_types
    contexts = {}
    bindings = {}
    preset = calculate_preset(ocel.graph.eog)

    log["event_objects"] = log.apply(lambda x: [(ot, o) for ot in object_types for o in x[ot]], axis=1)
    exploded_log = log.explode("event_objects")

    # Erstelle Argumentliste für den Pool
    args = [(event, preset[event], log, exploded_log, object_types) for event in preset.keys()]

    # Bestimme optimale Anzahl der Prozesse
    num_processes = min(len(preset), math.ceil(len(preset) / 1000))

    # Verwende Multiprocessing Pool
    with Pool(processes=5) as pool:
        results = list(tqdm(pool.imap(process_event, args), total=len(args)))

    # Verarbeite die Ergebnisse
    for event, context, binding_sequence in results:
        contexts[event] = context
        bindings[event] = binding_sequence

    return contexts, bindings
