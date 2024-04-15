#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from app import create_app
from flask_script import Manager, Shell

app = create_app('production')
manager = Manager(app=app)

def make_shell_context():
	return dict(app=app)

manager.add_command("shell", Shell(make_context=make_shell_context))


if __name__ == '__main__':
	manager.run()
