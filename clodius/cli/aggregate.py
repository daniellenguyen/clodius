# -*- coding: utf-8 -*-
from __future__ import division, print_function

from . import cli

import click
import clodius.tiles as ct
import collections as col
import h5py
import math
import negspy.coordinates as nc
import numpy as np
import os
import os.path as op
import pyBigWig as pbw
import random
import slugid
import sqlite3
import sys
import time

@cli.group()
def aggregate():
    '''
    Aggregate a data file so that it stores the data at multiple
    resolutions.
    '''
    pass

def store_meta_data(cursor, zoom_step, max_length, assembly, chrom_names, 
        chrom_sizes, tile_size, max_zoom, max_width):
    print("chrom_names:", chrom_names)

    cursor.execute('''
        CREATE TABLE tileset_info
        (
            zoom_step INT,
            max_length INT,
            assembly text,
            chrom_names text,
            chrom_sizes text,
            tile_size REAL,
            max_zoom INT,
            max_width REAL
        )
        ''')

    cursor.execute('INSERT INTO tileset_info VALUES (?,?,?,?,?,?,?,?)',
            (zoom_step, max_length, assembly, 
                "\t".join(chrom_names), "\t".join(map(str,chrom_sizes)),
                tile_size, max_zoom, max_width))
    cursor.commit()

    pass

# all entries are broked up into ((tile_pos), [entry]) tuples
# we just need to reduce the tiles so that no tile contains more than
# max_entries_per_tile entries
# (notice that [entry] is an array), this format will be important when
# reducing to the most important values
def reduce_values_by_importance(entry1, entry2, max_entries_per_tile=100, reverse_importance=False):
    def extract_key(entries):
        return [(e[-2], e) for e in entries]
    by_uid = dict(extract_key(entry1) + extract_key(entry2))
    combined_by_uid = by_uid.values()

    if reverse_importance:
        combined_entries = sorted(combined_by_uid,
                key=lambda x: float(x[-1]))
    else:
        combined_entries = sorted(combined_by_uid,
                key=lambda x: -float(x[-1]))

    byKey = {}

    return combined_entries[:max_entries_per_tile]

def _bedpe(filepath, output_file, assembly, importance_column, has_header, max_per_tile, 
        tile_size, max_zoom=None, chromosome=None, 
        chr1_col=0, from1_col=1, to1_col=2,
        chr2_col=3, from2_col=4, to2_col=5):
    print('output_file:', output_file)

    if filepath.endswith('.gz'):
        print("gzip")
        f = gzip.open(filepath, 'rt')
    else:
        print("plain")
        f = open(filepath, 'r')

    if output_file is None:
        output_file = filepath + ".multires.db"
    else:
        output_file = output_file

    if op.exists(output_file):
        os.remove(output_file)

    def line_to_dict(line):
        parts = line.split()
        d = {}
        try:
            d['xs'] = [nc.chr_pos_to_genome_pos(parts[chr1_col], int(parts[from1_col]), assembly), 
                          nc.chr_pos_to_genome_pos(parts[chr1_col], int(parts[to1_col]), assembly)]
            d['ys'] = [nc.chr_pos_to_genome_pos(parts[chr2_col], int(parts[from2_col]), assembly), 
                        nc.chr_pos_to_genome_pos(parts[chr2_col], int(parts[to2_col]), assembly)]
        except KeyError:
            error_str = ("ERROR converting chromosome position to genome position. "
                        "Please make sure you've specified the correct assembly "
                        "using the --assembly option. "
                        "Current assembly: {}, chromosomes: {},{}".format(assembly,
                    parts[chr1_col], parts[chr2_col]))
            raise(KeyError(error_str))

        d['uid'] = slugid.nice().decode('utf-8')

        d['chrOffset'] = d['xs'][0] - int(parts[from1_col])

        if importance_column is None:
            d['importance'] = max(d['xs'][1] - d['xs'][0], d['ys'][1] - d['ys'][0]) 
        elif importance_column == 'random':
            d['importance'] = random.random()
        else:
            d['importance'] = float(d[importance_column])

        d['fields'] = line

        return d

    entries = []

    if has_header:
        f.readline()
    else:
        first_line = f.readline().strip()
        try:
            parts = first_line.split()

            '''
            print("chr1_col", chr1_col, "chr2_col", chr2_col, 
                  "from1_col:", from1_col, "from2_col", from2_col, 
                  "to1_col", to1_col, "to2_col", to2_col)
            '''

            pos = int(parts[from1_col])
            pos = int(parts[to1_col])
            pos = int(parts[from2_col])
            pos = int(parts[to2_col])
        except ValueError as ve:
            error_str = "Couldn't convert one of the bedpe coordinates to an integer. If the input file contains a header, make sure to indicate that with the --has-header option. Line: {}".format(first_line)
            raise(ValueError(error_str))
        entries = [line_to_dict(first_line)]

    entries += [line_to_dict(line.strip()) for line in f]

    # We neeed chromosome information as well as the assembly size to properly
    # tile this data
    tile_size = tile_size
    chrom_info = nc.get_chrominfo(assembly)
    assembly_size = chrom_info.total_length+1
    #max_zoom = int(math.ceil(math.log(assembly_size / min_feature_width) / math.log(2)))
    max_zoom = int(math.ceil(math.log(assembly_size / tile_size) / math.log(2)))
    '''
    if max_zoom is not None and max_zoom < max_zoom:
        max_zoom = max_zoom
    '''

    # this script stores data in a sqlite database
    sqlite3.register_adapter(np.int64, lambda val: int(val))
    conn = sqlite3.connect(output_file)

    # store some meta data
    store_meta_data(conn, 1, 
            max_length = assembly_size,
            assembly = assembly,
            chrom_names = nc.get_chromorder(assembly),
            chrom_sizes = nc.get_chromsizes(assembly),
            tile_size = tile_size,
            max_zoom = max_zoom,
            max_width = tile_size * 2 ** max_zoom)

    max_width = tile_size * 2 ** max_zoom
    uid_to_entry = {}

    c = conn.cursor()
    c.execute(
    '''
    CREATE TABLE intervals
    (
        id int PRIMARY KEY,
        zoomLevel int,
        importance real,
        fromX int,
        toX int,
        fromY int,
        toY int,
        chrOffset int,
        uid text,
        fields text
    )
    ''')

    print("creating rtree")
    c.execute('''
        CREATE VIRTUAL TABLE position_index USING rtree(
            id,
            rFromX, rToX,
            rFromY, rToY
        )
        ''')

    curr_zoom = 0
    counter = 0
    
    max_viewable_zoom = max_zoom

    if max_zoom is not None and max_zoom < max_zoom:
        max_viewable_zoom = max_zoom

    tile_counts = col.defaultdict(lambda: col.defaultdict(lambda: col.defaultdict(int)))
    entries = sorted(entries, key=lambda x: -x['importance'])
    
    counter = 0
    for d in entries:
        curr_zoom = 0

        while curr_zoom <= max_zoom:
            tile_width = tile_size * 2 ** (max_zoom - curr_zoom)
            #print("d:", d)
            tile_from = list(map(lambda x: x / tile_width, [d['xs'][0], d['ys'][0]] ))
            tile_to = list(map(lambda x: x / tile_width, [d['xs'][1], d['ys'][1]]))

            empty_tiles = True

            # go through and check if any of the tiles at this zoom level are full

            for i in range(int(tile_from[0]), int(tile_to[0])+1):
                if not empty_tiles:
                    break

                for j in range(int(tile_from[1]), int(tile_to[1])+1):
                    if tile_counts[curr_zoom][i][j] > max_per_tile:

                        empty_tiles = False
                        break

            
            if empty_tiles:
                # they're all empty so add this interval to this zoom level
                for i in range(int(tile_from[0]), int(tile_to[0])+1):
                    for j in range(int(tile_from[1]), int(tile_to[1])+1):
                        tile_counts[curr_zoom][i][j] += 1

                #print("adding:", curr_zoom, d)
                exec_statement = 'INSERT INTO intervals VALUES (?,?,?,?,?,?,?,?,?,?)'
                ret = c.execute(
                        exec_statement,
                        (counter, curr_zoom, 
                            d['importance'],
                            d['xs'][0], d['xs'][1],
                            d['ys'][0], d['ys'][1],
                            d['chrOffset'], 
                            d['uid'],
                            d['fields'])
                        )
                conn.commit()

                exec_statement = 'INSERT INTO position_index VALUES (?,?,?,?,?)'
                ret = c.execute(
                        exec_statement,
                        (counter, d['xs'][0], d['xs'][1], 
                            d['ys'][0], d['ys'][1])  #add counter as a primary key
                        )
                conn.commit()

                counter += 1
                break

            curr_zoom += 1

    return

def _bedfile(filepath, output_file, assembly, importance_column, has_header, 
        chromosome, max_per_tile, tile_size, delimiter, chromsizes_filename,
        offset):
    if output_file is None:
        output_file = filepath + ".multires"
    else:
        output_file = output_file

    if op.exists(output_file):
        os.remove(output_file)

    bed_file = open(filepath, 'r')

    if chromsizes_filename is not None:
        chrom_info = nc.get_chrominfo_from_file(chromsizes_filename)
        chrom_names = chrom_info.chrom_order
        chrom_sizes = [chrom_info.chrom_lengths[c] for c in chrom_info.chrom_order]
    else:
        chrom_info = nc.get_chrominfo(assembly)
        chrom_names = nc.get_chromorder(assembly)
        chrom_sizes = nc.get_chromsizes(assembly)

    print("chrom_names:", chrom_info.chrom_order)
    print("chrom_sizes:", chrom_sizes)


    def line_to_np_array(line):
        '''
        Convert a bed file line to a numpy array which can later
        be used as an entry in an h5py file.
        '''
        try:
            start = int(line[1])
            stop = int(line[2])
        except ValueError:
            raise ValueError("Error parsing the position, line: {}".format(line))

        chrom = line[0]

        if importance_column is None:
            importance = stop - start
        elif importance_column == 'random':
            importance = random.random()
        else:
            importance = int(line[int(importance_column)-1])

        # convert chromosome coordinates to genome coordinates

        genome_start = chrom_info.cum_chrom_lengths[chrom] + start + offset
        #nc.chr_pos_to_genome_pos(str(chrom), start, assembly)
        genome_end = chrom_info.cum_chrom_lengths[chrom] + start + offset
        #nc.chr_pos_to_genome_pos(chrom, stop, assembly)

        pos_offset = genome_start - start
        parts = {
                    'startPos': genome_start,
                    'endPos': genome_end,
                    'uid': slugid.nice().decode('utf-8'),
                    'chrOffset': pos_offset,
                    'fields': '\t'.join(line),
                    'importance': importance,
                    'chromosome': str(chrom)
                    }

        return parts

    dset = []

    if has_header:
        bed_file.readline()
    else:
        line = bed_file.readline().strip()
        dset += [line_to_np_array(line.strip().split(delimiter))]

    for line in bed_file:
        dset += [line_to_np_array(line.strip().split(delimiter))]
    
    if chromosome is not None:
        dset = [d for d in dset if d['chromosome'] == chromosome]

    # We neeed chromosome information as well as the assembly size to properly
    # tile this data
    tile_size = tile_size

    #if chromosome is None:
    assembly_size = chrom_info.total_length+1
    '''
    else:
        try:
            assembly_size = chrom_info.chrom_lengths[chromosome]
        except KeyError:
            print("ERROR: Chromosome {} not found in assembly {}.".format(chromosome, assembly), file=sys.stderr)
            return 1
    '''

    #max_zoom = int(math.ceil(math.log(assembly_size / min_feature_width) / math.log(2)))
    max_zoom = int(math.ceil(math.log(assembly_size / tile_size) / math.log(2)))
    '''
    if max_zoom is not None and max_zoom < max_zoom:
        max_zoom = max_zoom
    '''

    # this script stores data in a sqlite database
    import sqlite3
    sqlite3.register_adapter(np.int64, lambda val: int(val))
    conn = sqlite3.connect(output_file)

    # store some meta data
    store_meta_data(conn, 1,
            max_length = assembly_size,
            assembly = assembly,
            chrom_names = chrom_names,
            chrom_sizes = chrom_sizes,
            tile_size = tile_size,
            max_zoom = max_zoom,
            max_width = tile_size * 2 ** max_zoom)

    max_width = tile_size * 2 ** max_zoom
    uid_to_entry = {}

    intervals = []

    # store each bed file entry as an interval
    for d in dset:
        uid = d['uid']
        uid_to_entry[uid] = d
        intervals += [(d['startPos'], d['endPos'], uid)]

    tile_width = tile_size

    removed = set()

    c = conn.cursor()
    c.execute(
    '''
    CREATE TABLE intervals
    (
        id int PRIMARY KEY,
        zoomLevel int,
        importance real,
        startPos int,
        endPos int,
        chrOffset int,
        uid text,
        fields text
    )
    ''')

    c.execute('''
        CREATE VIRTUAL TABLE position_index USING rtree(
            id,
            rStartPos, rEndPos
        )
        ''')

    curr_zoom = 0
    counter = 0

    max_viewable_zoom = max_zoom

    if max_zoom is not None and max_zoom < max_zoom:
        max_viewable_zoom = max_zoom

    while curr_zoom <= max_viewable_zoom and len(intervals) > 0:
        # at each zoom level, add the top genes
        tile_width = tile_size * 2 ** (max_zoom - curr_zoom)

        for tile_num in range(max_width // tile_width):
            # go over each tile and distribute the remaining values
            #values = interval_tree[tile_num * tile_width: (tile_num+1) * tile_width]
            from_value = tile_num * tile_width
            to_value = (tile_num + 1) * tile_width
            entries = [i for i in intervals if (i[0] < to_value and i[1] > from_value)]
            values_in_tile = sorted(entries,
                    key=lambda x: -uid_to_entry[x[-1]]['importance'])[:max_per_tile]   # the importance is always the last column
                                                            # take the negative because we want to prioritize
                                                            # higher values

            if len(values_in_tile) > 0:
                for v in values_in_tile:
                    counter += 1

                    value = uid_to_entry[v[-1]]

                    # one extra question mark for the primary key
                    exec_statement = 'INSERT INTO intervals VALUES (?,?,?,?,?,?,?,?)'
                    #print("value:", value['startPos'])

                    ret = c.execute(
                            exec_statement,
                            # primary key, zoomLevel, startPos, endPos, chrOffset, line
                            (counter, curr_zoom,
                                value['importance'],
                                value['startPos'], value['endPos'],
                                value['chrOffset'],
                                value['uid'],
                                value['fields'])
                            )
                    conn.commit()

                    exec_statement = 'INSERT INTO position_index VALUES (?,?,?)'
                    ret = c.execute(
                            exec_statement,
                            (counter, value['startPos'], value['endPos'])  #add counter as a primary key
                            )
                    conn.commit()
                    intervals.remove(v)
        #print ("curr_zoom:", curr_zoom, file=sys.stderr)
        curr_zoom += 1

    conn.commit()
    conn.close()

    return


def _bigwig(filepath, chunk_size=14, zoom_step=8, tile_size=1024, output_file=None, assembly='hg19', 
        chromsizes_filename=None, chromosome=None):
    last_end = 0
    data = []

    if output_file is None:
        if chromosome is None:
            output_file = op.splitext(filepath)[0] + '.hitile'
        else:
            output_file = op.splitext(filepath)[0] + '.' + chromosome + '.hitile'

    # Override the output file if it existts
    if op.exists(output_file):
        os.remove(output_file)
    f = h5py.File(output_file, 'w')

    if chromsizes_filename is not None:
        chrom_info = nc.get_chrominfo_from_file(chromsizes_filename)
        chrom_order = [a for a in nc.get_chromorder_from_file(chromsizes_filename)]
        chrom_sizes = nc.get_chromsizes_from_file(chromsizes_filename)
    else:
        print("there")
        chrom_info = nc.get_chrominfo(assembly)
        chrom_order = [a for a in nc.get_chromorder(assembly)]
        chrom_sizes = nc.get_chromsizes(assembly)

    print("chrom_order:", chrom_order)
    assembly_size = chrom_info.total_length

    tile_size = tile_size
    chunk_size = tile_size * 2**chunk_size     # how many values to read in at once while tiling

    dsets = []     # data sets at each zoom level
    nan_dsets = []

    # initialize the arrays which will store the values at each stored zoom level
    z = 0
    positions = []   # store where we are at the current dataset
    data_buffers = [[]]
    nan_data_buffers = [[]]

    while assembly_size / 2 ** z > tile_size:
        dset_length = math.ceil(assembly_size / 2 ** z)
        dsets += [f.create_dataset('values_' + str(z), (dset_length,), dtype='f',compression='gzip')]
        nan_dsets += [f.create_dataset('nan_values_' + str(z), (dset_length,), dtype='f',compression='gzip')]

        data_buffers += [[]]
        nan_data_buffers += [[]]

        positions += [0]
        z += zoom_step


    # load the bigWig file
    bwf = pbw.open(filepath)

    # store some meta data
    d = f.create_dataset('meta', (1,), dtype='f')

    if chromosome is not None:
        d.attrs['min-pos'] = chrom_info.cum_chrom_lengths[chromosome]
        d.attrs['max-pos'] = chrom_info.cum_chrom_lengths[chromosome] + bwf.chroms()[chromosome]
    else:
        d.attrs['min-pos'] = 0
        d.attrs['max-pos'] = assembly_size

    '''
    print("chroms.keys:", bwf.chroms().keys())
    print("chroms.values:", bwf.chroms().values())
    '''

    d.attrs['zoom-step'] = zoom_step
    d.attrs['max-length'] = assembly_size
    d.attrs['assembly'] = assembly
    d.attrs['chrom-names'] = [a.encode('utf-8') for a in chrom_order]
    d.attrs['chrom-sizes'] = chrom_sizes
    d.attrs['chrom-order'] = [a.encode('utf-8') for a in chrom_order]
    d.attrs['tile-size'] = tile_size
    d.attrs['max-zoom'] = max_zoom =  math.ceil(math.log(d.attrs['max-length'] / tile_size) / math.log(2))
    d.attrs['max-width'] = tile_size * 2 ** max_zoom
    d.attrs['max-position'] = 0

    print("assembly size (max-length)", d.attrs['max-length'])
    print("max-width", d.attrs['max-width'])
    print("max_zoom:", d.attrs['max-zoom'])
    print("chunk-size:", chunk_size)
    print("chrom-order", d.attrs['chrom-order'])

    t1 = time.time()

    curr_zoom = 0

    def add_values_to_data_buffers(buffers_to_add, nan_buffers_to_add):
        curr_zoom = 0

        data_buffers[0] += buffers_to_add
        nan_data_buffers[0] += nan_buffers_to_add

        curr_time = time.time() - t1
        percent_progress = (positions[curr_zoom] + 1) / float(assembly_size)
        print("position: {} progress: {:.2f} elapsed: {:.2f} remaining: {:.2f}".format(positions[curr_zoom] + 1, percent_progress,
            curr_time, curr_time / (percent_progress) - curr_time))

        while len(data_buffers[curr_zoom]) >= chunk_size:
            # get the current chunk and store it, converting nans to 0
            print("len(data_buffers[curr_zoom])", len(data_buffers[curr_zoom]))
            curr_chunk = np.array(data_buffers[curr_zoom][:chunk_size])
            nan_curr_chunk = np.array(nan_data_buffers[curr_zoom][:chunk_size])
            #curr_chunk[np.isnan(curr_chunk)] = 0
            '''
            print("1cc:", sum(curr_chunk))
            print("1db:", data_buffers[curr_zoom][:chunk_size])
            print("1curr_chunk:", nan_curr_chunk)
            '''
            print("positions[curr_zoom]:", positions[curr_zoom])

            dsets[curr_zoom][positions[curr_zoom]:positions[curr_zoom]+chunk_size] = curr_chunk
            nan_dsets[curr_zoom][positions[curr_zoom]:positions[curr_zoom]+chunk_size] = nan_curr_chunk

            # aggregate nan values
            #nan_curr_chunk[np.isnan(curr_chunk)] = 0
            #print("1na_cc:", sum(nan_curr_chunk))

            # aggregate and store aggregated values in the next zoom_level's data
            data_buffers[curr_zoom+1] += list(ct.aggregate(curr_chunk, 2 ** zoom_step))
            nan_data_buffers[curr_zoom+1] += list(ct.aggregate(nan_curr_chunk, 2 ** zoom_step))

            data_buffers[curr_zoom] = data_buffers[curr_zoom][chunk_size:]
            nan_data_buffers[curr_zoom] = nan_data_buffers[curr_zoom][chunk_size:]

            data = data_buffers[curr_zoom+1]
            nan_data = nan_data_buffers[curr_zoom+1]

            # do the same for the nan values buffers

            positions[curr_zoom] += chunk_size
            curr_zoom += 1

            if curr_zoom * zoom_step >= max_zoom:
                break

    # Do we only want values from a single chromosome?
    if chromosome is not None:
        chroms_to_use = [chromosome]
    else:
        chroms_to_use = chrom_order

    for chrom in chroms_to_use:
        print("chrom:", chrom)
        '''
        if chrom not in bwf.chroms():
            print("skipping chrom (not in bigWig file):",
            chrom, chrom_info.chrom_lengths[chrom])
            continue
        '''

        counter = 0
        # chrom_size = bwf.chroms()[chrom]
        chrom_size = chrom_info.chrom_lengths[chrom]

        # print("chrom_size:", chrom_size, bwf.chroms()[chrom])
        d.attrs['max-position'] += chrom_size

        while counter < chrom_size:
            remaining = min(chunk_size, chrom_size - counter)

            if chrom not in bwf.chroms():
                values = [np.nan] * remaining
                nan_values = [1] * remaining
            else:
                values = bwf.values(chrom, counter, counter + remaining)
                nan_values = np.isnan(values).astype('i4')

            # print("counter:", counter, "remaining:", remaining,
            # "counter + remaining:", counter + remaining)
            counter += remaining
            curr_zoom = 0

            add_values_to_data_buffers(list(values), list(nan_values))

    while True:
        # get the current chunk and store it
        chunk_size = len(data_buffers[curr_zoom])
        curr_chunk = np.array(data_buffers[curr_zoom][:chunk_size])
        nan_curr_chunk = np.array(nan_data_buffers[curr_zoom][:chunk_size])

        dsets[curr_zoom][positions[curr_zoom]:positions[curr_zoom]+chunk_size] = curr_chunk
        nan_dsets[curr_zoom][positions[curr_zoom]:positions[curr_zoom]+chunk_size] = nan_curr_chunk

        # aggregate and store aggregated values in the next zoom_level's data
        data_buffers[curr_zoom+1] += list(ct.aggregate(curr_chunk, 2 ** zoom_step))
        nan_data_buffers[curr_zoom+1] += list(ct.aggregate(nan_curr_chunk, 2 ** zoom_step))

        data_buffers[curr_zoom] = data_buffers[curr_zoom][chunk_size:]
        nan_data_buffers[curr_zoom] = nan_data_buffers[curr_zoom][chunk_size:]

        data = data_buffers[curr_zoom+1]
        nan_data = nan_data_buffers[curr_zoom+1]

        positions[curr_zoom] += chunk_size
        curr_zoom += 1

        # we've created enough tile levels to cover the entire maximum width
        if curr_zoom * zoom_step >= max_zoom:
            break

    # still need to take care of the last chunk

    data = np.array(data)
    t1 = time.time()
    pass

##################################################################################################
def _bedgraph(filepath, output_file, assembly, chrom_col, 
        from_pos_col, to_pos_col, value_col, has_header, 
        chromosome, tile_size, chunk_size, method, nan_value,
        transform, count_nan, chromsizes_filename, zoom_step):
    last_end = 0
    data = []

    if output_file is None:
        output_file = op.splitext(filepath)[0] + '.hitile'

    print("output file:", output_file)

    # Override the output file if it existts
    if op.exists(output_file):
        os.remove(output_file)
    f = h5py.File(output_file, 'w')

    # get the information about the chromosomes in this assembly
    if chromsizes_filename is not None:
        chrom_info = nc.get_chrominfo_from_file(chromsizes_filename)
        chrom_order = [a.encode('utf-8') for a in nc.get_chromorder_from_file(chromsizes_filename)]
        chrom_sizes = nc.get_chromsizes_from_file(chromsizes_filename)
    else:
        chrom_info = nc.get_chrominfo(assembly)
        chrom_order = [a.encode('utf-8') for a in nc.get_chromorder(assembly)]
        chrom_sizes = nc.get_chromsizes(assembly)

    assembly_size = chrom_info.total_length
    print('assembly_size:', assembly_size)

    tile_size = tile_size
    chunk_size = tile_size * 2**chunk_size     # how many values to read in at once while tiling

    dsets = []     # data sets at each zoom level
    nan_dsets = []  # store nan values

    # initialize the arrays which will store the values at each stored zoom level
    z = 0
    positions = []   # store where we are at the current dataset
    data_buffers = [[]]
    nan_data_buffers = [[]]

    while assembly_size / 2 ** z > tile_size:
        dset_length = math.ceil(assembly_size / 2 ** z)
        dsets += [f.create_dataset('values_' + str(z), (dset_length,), dtype='f',compression='gzip')]
        nan_dsets += [f.create_dataset('nan_values_' + str(z), (dset_length,), dtype='f',compression='gzip')]

        data_buffers += [[]]
        nan_data_buffers += [[]]

        positions += [0]
        z += zoom_step

    #print("dsets[0][-10:]", dsets[0][-10:])

    # load the bigWig file
    #print("filepath:", filepath)

    # store some meta data
    d = f.create_dataset('meta', (1,), dtype='f')

    print("assembly:", assembly)
    #print("chrom_info:", nc.get_chromorder(assembly))

    d.attrs['zoom-step'] = zoom_step
    d.attrs['max-length'] = assembly_size
    d.attrs['assembly'] = assembly
    d.attrs['chrom-names'] = chrom_order
    d.attrs['chrom-sizes'] = chrom_sizes
    d.attrs['chrom-order'] = chrom_order
    d.attrs['tile-size'] = tile_size
    d.attrs['max-zoom'] = max_zoom =  math.ceil(math.log(d.attrs['max-length'] / tile_size) / math.log(2))
    d.attrs['max-width'] = tile_size * 2 ** max_zoom
    d.attrs['max-position'] = 0

    print("assembly size (max-length)", d.attrs['max-length'])
    print("max-width", d.attrs['max-width'])
    print("max_zoom:", d.attrs['max-zoom'])
    print("chunk-size:", chunk_size)
    print("chrom-order", d.attrs['chrom-order'])

    t1 = time.time()

    # are we reading the input from stdin or from a file?

    if filepath == '-':
        f = sys.stdin
    else:
        if filepath.endswith('.gz'):
            import gzip
            f = gzip.open(filepath, 'rt')
        else:
            f = open(filepath, 'r')

    curr_zoom = 0

    def add_values_to_data_buffers(buffers_to_add, nan_buffers_to_add):
        curr_zoom = 0

        data_buffers[0] += buffers_to_add
        nan_data_buffers[0] += nan_buffers_to_add

        curr_time = time.time() - t1
        percent_progress = (positions[curr_zoom] + 1) / float(assembly_size)
        print("position: {} progress: {:.2f} elapsed: {:.2f} remaining: {:.2f}".format(positions[curr_zoom] + 1, percent_progress,
            curr_time, curr_time / (percent_progress) - curr_time))

        while len(data_buffers[curr_zoom]) >= chunk_size:
            # get the current chunk and store it, converting nans to 0
            print("len(data_buffers[curr_zoom])", len(data_buffers[curr_zoom]))
            curr_chunk = np.array(data_buffers[curr_zoom][:chunk_size])
            nan_curr_chunk = np.array(nan_data_buffers[curr_zoom][:chunk_size])
            #curr_chunk[np.isnan(curr_chunk)] = 0
            '''
            print("1cc:", sum(curr_chunk))
            print("1db:", data_buffers[curr_zoom][:chunk_size])
            print("1curr_chunk:", nan_curr_chunk)
            '''
            print("positions[curr_zoom]:", positions[curr_zoom])

            dsets[curr_zoom][positions[curr_zoom]:positions[curr_zoom]+chunk_size] = curr_chunk
            nan_dsets[curr_zoom][positions[curr_zoom]:positions[curr_zoom]+chunk_size] = nan_curr_chunk

            # aggregate nan values
            #nan_curr_chunk[np.isnan(curr_chunk)] = 0
            #print("1na_cc:", sum(nan_curr_chunk))

            # aggregate and store aggregated values in the next zoom_level's data
            data_buffers[curr_zoom+1] += list(ct.aggregate(curr_chunk, 2 ** zoom_step))
            nan_data_buffers[curr_zoom+1] += list(ct.aggregate(nan_curr_chunk, 2 ** zoom_step))

            data_buffers[curr_zoom] = data_buffers[curr_zoom][chunk_size:]
            nan_data_buffers[curr_zoom] = nan_data_buffers[curr_zoom][chunk_size:]

            data = data_buffers[curr_zoom+1]
            nan_data = nan_data_buffers[curr_zoom+1]

            # do the same for the nan values buffers

            positions[curr_zoom] += chunk_size
            curr_zoom += 1

            if curr_zoom * zoom_step >= max_zoom:
                break


    values = []
    nan_values = []

    if has_header:
        f.readline()

    # the genome position up to which we've filled in values
    curr_genome_pos = 0

    # keep track of the previous value so that we can use it to fill in NAN values
    prev_value = 0

    for line in f:
        # each line should indicate a chromsome, start position and end position
        parts = line.strip().split()

        start_genome_pos = chrom_info.cum_chrom_lengths[parts[chrom_col-1]] + int(parts[from_pos_col-1])         
        #print("len(values):", len(values), curr_genome_pos, start_genome_pos)
        #print("line:", line)

        if start_genome_pos - curr_genome_pos > 1:
            values += [np.nan] * (start_genome_pos - curr_genome_pos - 1)
            nan_values += [1] * (start_genome_pos - curr_genome_pos - 1)

            curr_genome_pos += (start_genome_pos - curr_genome_pos - 1)


        # count how many nan values there are in the dataset
        nan_count = 1 if parts[value_col-1] == nan_value else 0

        # if the provided values are log2 transformed, we have to un-transform them
        if transform == 'exp2':
            value = 2 ** float(parts[value_col-1]) if not parts[value_col-1] == nan_value else np.nan
        else:
            value = float(parts[value_col-1]) if not parts[value_col-1] == nan_value else np.nan


        # print("pos:", int(parts[to_pos_col-1]) - int(parts[from_pos_col-1]))
        # we're going to add as many values are as specified in the bedfile line
        values_to_add = [value] * (int(parts[to_pos_col-1]) - int(parts[from_pos_col-1]))
        nan_counts_to_add = [nan_count] * (int(parts[to_pos_col-1]) - int(parts[from_pos_col-1]))
        
        values += values_to_add
        nan_values += nan_counts_to_add

        d.attrs['max-position'] = start_genome_pos + len(values_to_add) 

        #print("values:", values[:30])

        curr_genome_pos += len(values_to_add)

        while len(values) > chunk_size:
            print("len(values):", len(values), chunk_size)
            print("line:", line)
            add_values_to_data_buffers(values[:chunk_size], nan_values[:chunk_size])
            values = values[chunk_size:]
            nan_values = nan_values[chunk_size:]


    add_values_to_data_buffers(values, nan_values)

    # store the remaining data
    while True:
        # get the current chunk and store it
        chunk_size = len(data_buffers[curr_zoom])
        curr_chunk = np.array(data_buffers[curr_zoom][:chunk_size])
        nan_curr_chunk = np.array(nan_data_buffers[curr_zoom][:chunk_size])

        '''
        print("2curr_chunk", curr_chunk)
        print("2curr_zoom:", curr_zoom)
        print("2db", data_buffers[curr_zoom][:100])
        '''

        dsets[curr_zoom][positions[curr_zoom]:positions[curr_zoom]+chunk_size] = curr_chunk
        nan_dsets[curr_zoom][positions[curr_zoom]:positions[curr_zoom]+chunk_size] = nan_curr_chunk

        #print("chunk_size:", chunk_size, "len(curr_chunk):", len(curr_chunk), "len(nan_curr_chunk)", len(nan_curr_chunk))

        # aggregate and store aggregated values in the next zoom_level's data
        data_buffers[curr_zoom+1] += list(ct.aggregate(curr_chunk, 2 ** zoom_step))
        nan_data_buffers[curr_zoom+1] += list(ct.aggregate(nan_curr_chunk, 2 ** zoom_step))

        data_buffers[curr_zoom] = data_buffers[curr_zoom][chunk_size:]
        nan_data_buffers[curr_zoom] = nan_data_buffers[curr_zoom][chunk_size:]

        data = data_buffers[curr_zoom+1]
        nan_data = nan_data_buffers[curr_zoom+1]

        positions[curr_zoom] += chunk_size
        curr_zoom += 1

        # we've created enough tile levels to cover the entire maximum width
        if curr_zoom * zoom_step >= max_zoom:
            break

    # still need to take care of the last chunk

@aggregate.command()
@click.argument(
        'filepath',
        metavar='FILEPATH'
        )
@click.option(
        '--output-file',
        '-o',
        default=None,
        help="The default output file name to use. If this isn't"
             "specified, clodius will replace the current extension"
             "with .hitile"
        )
@click.option(
        '--assembly',
        '-a',
        help='The genome assembly that this file was created against',
        type=click.Choice(nc.available_chromsizes()),
        default='hg19')
@click.option(
        '--chromosome',
        default=None,
        help="Only extract values for a particular chromosome."
             "Use all chromosomes if not set.")
@click.option(
        '--tile-size',
        '-t',
        default=1024,
        help="The number of data points in each tile."
             "Used to determine the number of zoom levels"
             "to create.")
@click.option(
        '--chunk-size',
        '-c',
        help='How many values to aggregate at once.'
             'Specified as a power of two multiplier of the tile'
             'size',
        default=14)
@click.option(
        '--chromosome-col',
        help="The column number (1-based) which contains the chromosome "
              "name",
        default=1)
@click.option(
        '--from-pos-col',
        help="The column number (1-based) which contains the starting "
             "position",
        default=2)
@click.option(
        '--to-pos-col',
        help="The column number (1-based) which contains the ending"
             "position",
        default=3)
@click.option(
        '--value-col',
        help="The column number (1-based) which contains the actual value"
             "position",
        default=4)
@click.option(
        '--has-header/--no-header',
        help="Does this file have a header that we should ignore",
        default=False)
@click.option(
        '--method',
        help='The method used to aggregate values (e.g. sum, average...)',
        type=click.Choice(['sum', 'average']),
        default='sum')
@click.option(
        '--nan-value',
        help='The string to use as a NaN value',
        type=str,
        default=None)
@click.option(
        '--transform',
        help='The method used to aggregate values (e.g. sum, average...)',
        type=click.Choice(['none', 'exp2']),
        default='none')
@click.option(
        '--count-nan',
        help="Simply count the number of nan values in the file",
        is_flag=True)
@click.option(
        '--chromsizes-filename',
        help="A file containing chromosome sizes and order",
        default=None)
@click.option(
        '--zoom-step',
        '-z',
        help="The number of intermediate aggregation levels to"
             "omit",
        default=8)
def bedgraph(filepath, output_file, assembly, chromosome_col, 
        from_pos_col, to_pos_col, value_col, has_header, 
        chromosome, tile_size, chunk_size, method, nan_value, 
        transform, count_nan, chromsizes_filename, zoom_step):
    _bedgraph(filepath, output_file, assembly, chromosome_col, 
        from_pos_col, to_pos_col, value_col, has_header, 
        chromosome, tile_size, chunk_size, method, nan_value, 
        transform, count_nan, chromsizes_filename, zoom_step)

@aggregate.command()
@click.argument(
        'filepath',
        metavar='FILEPATH'
        )
@click.option(
        '--output-file',
        '-o',
        default=None,
        help="The default output file name to use. If this isn't"
             "specified, clodius will replace the current extension"
             "with .hitile"
        )
@click.option(
        '--assembly',
        '-a',
        help='The genome assembly that this file was created against',
        default='hg19')
@click.option(
        '--chromosome',
        default=None,
        help="Only extract values for a particular chromosome."
             "Use all chromosomes if not set.")
@click.option(
        '--tile-size',
        '-t',
        default=1024,
        help="The number of data points in each tile."
             "Used to determine the number of zoom levels"
             "to create.")
@click.option(
        '--chunk-size',
        '-c',
        help='How many values to aggregate at once.'
             'Specified as a power of two multiplier of the tile'
             'size',
        default=14)
@click.option(
        '--chromsizes-filename',
        help="A file containing chromosome sizes and order",
        default=None)
@click.option(
        '--zoom-step',
        '-z',
        help="The number of intermediate aggregation levels to"
             "omit",
        default=8)
def bigwig(filepath, output_file, assembly, chromosome, tile_size, chunk_size, chromsizes_filename, zoom_step):
    _bigwig(filepath, chunk_size, zoom_step, tile_size, output_file, assembly, chromsizes_filename, chromosome)

@aggregate.command()
@click.argument( 
        'filepath',
        metavar='FILEPATH')
@click.option(
        '--output-file',
        '-o',
        default=None,
        help="The default output file name to use. If this isn't"
             "specified, clodius will replace the current extension"
             "with .multires.bed"
        )
@click.option(
        '--assembly',
        '-a',
        help='The genome assembly that this file was created against',
        default='hg19')
@click.option(
        '--importance-column',
        help='The column (1-based) containing information about how important'
        "that row is. If it's absent, then use the length of the region."
        "If the value is equal to `random`, then a random value will be"
        "used for the importance (effectively leading to random sampling)"
        )
@click.option(
        '--has-header/--no-header',
        help="Does this file have a header that we should ignore",
        default=False)
@click.option(
        '--chromosome',
        default=None,
        help="Only extract values for a particular chromosome."
             "Use all chromosomes if not set."
             )
@click.option(
        '--max-per-tile',
        default=100,
        type=int)
@click.option(
        '--tile-size', 
        default=1024,
        help="The number of nucleotides that the highest resolution tiles should span."
             "This determines the maximum zoom level"
        )
@click.option(
        '--delimiter',
        default=None,
        type=str)
@click.option(
        '--chromsizes-filename',
        help="A file containing chromosome sizes and order",
        default=None)
@click.option(
        '--offset',
        help="Apply an offset to all the coordinates in this file",
        type=int,
        default=0)
def bedfile(filepath, output_file, assembly, importance_column, has_header, 
        chromosome, max_per_tile, tile_size, delimiter, chromsizes_filename,
        offset):
    _bedfile(filepath, output_file, assembly, importance_column, has_header, 
            chromosome, max_per_tile, tile_size, delimiter, chromsizes_filename,
            offset)

@aggregate.command()
@click.argument( 
        'filepath',
        metavar='FILEPATH')
@click.option(
        '--output-file',
        '-o',
        default=None,
        help="The default output file name to use. If this isn't"
             "specified, clodius will replace the current extension"
             "with .bed2db"
        )
@click.option(
        '--assembly',
        '-a',
        help='The genome assembly that this file was created against',
        default='hg19')
@click.option(
        '--importance-column',
        help='The column (1-based) containing information about how important'
        "that row is. If it's absent, then use the length of the region."
        "If the value is equal to `random`, then a random value will be"
        "used for the importance (effectively leading to random sampling)",
        default='random'
        )
@click.option(
        '--has-header/--no-header',
        help="Does this file have a header that we should ignore",
        default=False)
@click.option(
        '--max-per-tile',
        default=100,
        type=int)
@click.option(
        '--tile-size', 
        default=1024,
        help="The number of nucleotides that the highest resolution tiles should span."
             "This determines the maximum zoom level"
        )
@click.option(
        '--chromosome',
        default=None,
        help="Only extract values for a particular chromosome."
             "Use all chromosomes if not set."
             )
@click.option(
        '--chr1-col',
        default=1,
        help="The column containing the first chromosome"
             )
@click.option(
        '--chr2-col',
        default=4,
        help="The column containing the second chromosome"
             )
@click.option(
        '--from1-col',
        default=2,
        help="The column containing the first start position"
             )
@click.option(
        '--from2-col',
        default=5,
        help="The column containing the second start position"
             )
@click.option(
        '--to1-col',
        default=3,
        help="The column containing the first end position"
             )
@click.option(
        '--to2-col',
        default=6,
        help="The column containing the second end position"
             )

def bedpe(filepath, output_file, assembly, importance_column, 
        has_header, max_per_tile, tile_size, chromosome,
        chr1_col, from1_col, to1_col,
        chr2_col, from2_col, to2_col):

    print("## chr1_col", chr1_col, "chr2_col", chr2_col, 
          "from1_col:", from1_col, "from2_col", from2_col, 
          "to1_col", to1_col, "to2_col", to2_col)
    _bedpe(filepath, output_file, assembly, importance_column, has_header, 
            max_per_tile, tile_size, chromosome,
            chr1_col=chr1_col-1, from1_col=from1_col-1, to1_col=to1_col-1,
            chr2_col=chr2_col-1, from2_col=from2_col-1, to2_col=to2_col-1
            )
