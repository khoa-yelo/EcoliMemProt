def find_values_by_key(dictionary, target_key):
    values = []
    for key, value in dictionary.items():
        if key == target_key:
            values.append(value)
        if isinstance(value, dict):
            values.extend(find_values_by_key(value, target_key))
    return values


def join_list(input_list):
    result = []
    for element in input_list:
        if isinstance(element, list):
            result.extend(element)
        else:
            result.append(element)
    return result
