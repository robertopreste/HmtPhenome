#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste


class Node:
    """
    This class is used to create nodes that are used to draw the network.
    """

    def __init__(self, entity: str, nid: int, label: str):
        self.entity = entity
        self.nid = nid  # node ID (1, 2, 3, etc)
        self.label = label
        self.uid = hash(self.entity + self.label)
