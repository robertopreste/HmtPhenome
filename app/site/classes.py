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
        self._background = None
        self._border = None
        self._connected = None

    @property
    def background(self):
        if self.entity == "v":
            self._background = "#F9CF45"
        elif self.entity == "g":
            self._background = "#739E82"
        elif self.entity == "d":
            self._background = "#D7816A"
        elif self.entity == "p":
            self._background = "#93B5C6"
        return self._background

    @property
    def border(self):
        if self.entity == "v":
            self._border = "#F9CF45"
        elif self.entity == "g":
            self._border = "#739E82"
        elif self.entity == "d":
            self._border = "#D7816A"
        elif self.entity == "p":
            self._border = "#93B5C6"
        return self._border

