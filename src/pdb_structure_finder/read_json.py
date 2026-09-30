import json


def json_data(json_file_path):
    try:
        with open(json_file_path, "r", encoding="utf-8") as file:
            content = json.load(file)
            return content
    except FileNotFoundError:
        raise
    except json.JSONDecodeError:
        json_structure()
        raise


# TODO insert a function for printing how a json shuold be formatted and suggest the use of the tool for the interactive construction of json
# NOTE Could be nice if the function would be called by the raise of an error in the reading file
def json_structure():

    return
