# -*- coding:utf-8 -*-


from flask import Blueprint


network = Blueprint('network', __name__)
from . import views