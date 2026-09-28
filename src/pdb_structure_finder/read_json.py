import json


def json_data(json_file_path):
    with open(json_file_path, "r", encoding="utf-8") as file:
        return json.load(json_file_path)


# FIXME: need to fix json.load error
