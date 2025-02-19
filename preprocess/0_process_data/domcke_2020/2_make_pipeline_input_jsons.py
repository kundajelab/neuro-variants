import os
import json

with open("pipeline.input.template.json", "r") as template_file:
    template_json = json.load(template_file)

tagalign_dir = '/oak/stanford/groups/akundaje/projects/neuro-variants/data/processed/domcke_2020/tagaligns'
output_dir = 'pipeline_input_jsons'

os.makedirs(output_dir, exist_ok=True)

for file_name in os.listdir(tagalign_dir):
    sample_name = '.'.join(file_name.split('.')[:3])
    print(sample_name)

    json_data = template_json.copy()
    json_data['atac.title'] = sample_name
    json_data["atac.tas"] = [os.path.join(json_data["atac.tas"][0], file_name)]

    output_file = os.path.join(output_dir, sample_name + '.input.json')
    with open(output_file, "w") as json_file:
        json.dump(json_data, json_file, indent=4)

