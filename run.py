#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from app import app


if __name__ == '__main__':
    app.run(debug=True, port=8080)  # remove this in production
