import os
import re
import shutil
import urllib.request
import urllib.error
import zipfile
#import rpy2.robjects.packages

import rpy2.robjects as robjects
import numpy as np
import networkx as nx

from bs4 import BeautifulSoup
from concurrent.futures import ThreadPoolExecutor
from functools import reduce


def get_movie_dataset(datadir):
    zip_url = 'https://dataverse.harvard.edu/api/access/datafile/:'\
              'persistentId?persistentId=doi:10.7910/DVN/T4HBA3/BIOFH6'
    tmp_zip, _ = urllib.request.urlretrieve(zip_url)

    if os.path.exists(datadir):
        shutil.rmtree(datadir)
    os.mkdir(datadir)

    # https://stackoverflow.com/a/4917469
    with zipfile.ZipFile(tmp_zip) as zip_file:
        for member in zip_file.namelist():
            filename = os.path.basename(member)
            # skip directories
            if not filename or filename[0] == '.':
                continue

            # copy file (taken from zipfile's extract)
            source = zip_file.open(member)
            target = open(os.path.join(datadir, filename), "wb")
            with source, target:
                shutil.copyfileobj(source, target)


def get_movie_info():
    MG_ordered_url = 'https://moviegalaxies.com/discover/movies/all/'\
                     '?order_by=title&ordering=asc&page='

    def MG_page_request(page_ind):
        req = urllib.request.Request(MG_ordered_url+str(page_ind), 
                                     headers={'User-Agent': 'Chrome Mozilla'})
        names_page = urllib.request.urlopen(req).read()
        soup = BeautifulSoup(names_page, 'html.parser')
        name_list = [el.text.split(os.linesep)[1] for el in soup.select('div[class=p-4]')]
        
        return name_list


    with ThreadPoolExecutor(max_workers=10) as executor:
        movie_info = executor.map(MG_page_request, range(1,27))
    movie_info = reduce(lambda x,y: x+y, movie_info)


    # Fixes names to be similar to MG Dataset
    str_fix = lambda s: s.replace(os.linesep, '').split('(')[0].strip(
                                  ' ').replace('.', '').upper()

    def IMDB_page_request(movie_name):
        # Get Movie ID
        search_page = urllib.request.urlopen("http://www.imdb.com/find?q=" + \
                                     f"{urllib.parse.quote(movie_name)}").read()
        soup = BeautifulSoup(search_page, 'html.parser')
        movie_id = soup.find('a', text=re.compile(movie_name))
        if movie_id is None or movie_id == []:
            return [[],[],[]]
        
        movie_id = movie_id.get('href')

        # Get Movie Genres
        movie_page = urllib.request.urlopen("https://www.imdb.com" + movie_id)
        soup = BeautifulSoup(movie_page, 'html.parser')
        genres = [el.text for el in soup.select('div[class=subtext] > '
                                                '[href^="/search/title?genres"]')]

        # Get Movies Characters
        try:
            chars_page = urllib.request.urlopen("https://www.imdb.com" + movie_id + "fullcredits")
            chars_page = chars_page.read()
            soup = BeautifulSoup(chars_page, 'html.parser')
            chars = [str_fix(el.text) for el in soup.select('.character')]
        except urllib.error.HTTPError:
            chars = []
        
        return [movie_name, genres, chars]


    with ThreadPoolExecutor(max_workers=32) as executor:
        movie_info = executor.map(IMDB_page_request, movie_info)


    return list(movie_info)


def rank_and_matching(sub_list, full_list):
    n1,n2 = len(sub_list), len(full_list)
    if n2 == 0:
        return [], [], False # no cast listed

    matching = [0]*n1
    rank = [0]*n1
    sub_opts = list(range(n1))

    # Check for main character
    k = 1
    for j in sub_opts:
        x = sub_list[j]
        if x == full_list[0] or x in full_list[0].split(' '):
            matching[j] = 1
            rank[j] = k

            k += 1
            sub_opts.remove(j)
            break # partial matches can keep matching

    if len(sub_opts) == n1:
        return [], [], False # negative result


    # Continue gathering ranks
    for i in range(1,n2):
        y = full_list[i]
        for j in sub_opts:
            x = sub_list[j]
            if x == y or x in y.split(' '):
                matching[j] = 1
                rank[j] = k

                k += 1
                sub_opts.remove(j)
                break # partial matches can keep matching

    return rank, matching, True


def match_and_convert(datadir, movie_info, min_matching):
    #rutils = rpy2.robjects.packages.importr('utils')
    #rutils.install_packages('network')
    #rnetwork = rpy2.robjects.packages.importr('network', on_conflict='warn')

    # Sort gexf_file by numerical value
    for gexf_file in sorted(os.listdir(datadir), key=lambda x: int(x.split('.')[0])):
        movie_net = nx.read_gexf(os.path.join(datadir, gexf_file))
        movie_nodes = movie_net.nodes

        char_names = [movie_nodes[el]['label'] for el in movie_nodes]
        n = len(movie_nodes)
        for i in range(len(movie_info)):
            all_chars = movie_info[i][2]
            if type(all_chars) != list:
                continue

            ranks, matching, has_main_char = rank_and_matching(char_names, all_chars)

            if has_main_char and sum(matching)/n > min_matching:
                rmat = robjects.r.matrix(robjects.IntVector(np.array(nx.adjacency_matrix(movie_net, 
                                         movie_nodes).todense()).flatten()), nrow=n, ncol=n)

                # Character ranks by importance
                k = max(ranks)
                for j in range(n):
                    if ranks[j] == 0:
                        ranks[j] = k + 1
                        k = k + 1
                
                #rnet = rnetwork.as_network(rmat, directed=False, matrix_type='a',
                #                           ignore_eval=False, names_eval='e_weights')
                #rnetwork.set_vertex_attribute(rnet, 'char_ranks', robjects.IntVector(ranks))
                #rnetwork.set_vertex_attribute(rnet, 'char_names', robjects.StrVector(char_names))

                movie_info[i][2] = rmat
                movie_info[i].append(robjects.IntVector(ranks))
                movie_info[i].append(robjects.StrVector(char_names))
                break

    return movie_info

import time
def save_rdata(datadir='data', min_matching=0.3):
    # Get MovieGalaxies Dataset
    get_movie_dataset(datadir)

    # Scrape Genres and Ordered (By Importance) Character List
    movie_info = get_movie_info()

    # Account for Discrepancies Between MG Dataset and Scraped Names
    # Also converts network to appropriate R object
    movie_info = match_and_convert(datadir, movie_info, min_matching)

    coll_dict = {}
    n = len(movie_info)
    for movie_list in movie_info:
        if type(movie_list[2]) == list:
                n -= 1
                continue

        for genre in movie_list[1]:
            if genre not in coll_dict:
                coll_dict[genre] = {}

            coll_dict[genre][movie_list[0]] = robjects.ListVector({'adj_mat':movie_list[2],
                                                'ranks':movie_list[3], 'names':movie_list[4]})

    coll_dict = robjects.ListVector({k:robjects.ListVector(v) for k,v in coll_dict.items()})
    robjects.r.assign('movies_by_genre', coll_dict)
    robjects.r("save(movies_by_genre, file='{}')".format('movie_data.RData'))
    print('Total Number Of Unique Movies:', n)

save_rdata()