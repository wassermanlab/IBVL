import os
from result import *
from pprint import pprint

from pyopencga.opencga_config import ClientConfiguration
from pyopencga.opencga_client import OpenCGAClient



def get_individual_details(resp):
    if not resp['annotationSets']: # check that this list is not empty
        return Err("Not found")      
    res = resp['annotationSets'][0]['annotations']
    return Ok(res)


if __name__=='__main__':
    # Configure and Connect
    config = ClientConfiguration({
            "rest": {
                    "host": os.environ.get("OPENCGA_HOST")
            }
    })
    oc = OpenCGAClient(config)

    # Enter your username and password here
    oc.login(os.environ.get("OPENCGA_USERNAME"), os.environ.get("OPENCGA_PASSWORD"))
    users = oc.users
    projects = oc.projects
    studies = oc.studies
    files = oc.files
    jobs = oc.jobs
    families = oc.families
    individuals = oc.individuals
    samples = oc.samples
    cohorts = oc.cohorts
    panels = oc.panels

    sample_result = oc.samples.search(study='1kG_phase3', limit=3, include='id')

    print(sample_result.responses[0]['result'])


    # See available methods

    #[res['id'] for res in sample_result.results()]
    for res in sample_result.results():
        print(res['id'])

    #[res for res in sample_result.results()]
    for res in sample_result.results():
        print(res)

    # Update Individual Metadata
    individuals.update(query_id="HG00099",data={"samples":["HG00099"],"sex":"MALE","karyotypicSex":"XY","lifeStatus":"ALIVE"})
    samples.update(query_id='HG00096',study='1kG_phase3',data={"individual":"HG00096","description":"demo description"} )
    res = samples.search(study="1kG_phase3",limit=5)
    #[(r['id'],) for r in res.results()]
    for r in res.results():
        print(r['id'])

    for r in res.results():
        print(r['id'], r.get('attributes').get('OPENCGA_INDIVIDUAL'))

    # VariablSet annotations
    # We can add custom annotations to individuals, samples, cohorts, etc through the idea of VariableSet 
    # similar to a table but more flexible
    individuals.update(query_id="HG00103",
                    data={"sex":"MALE","karyotypicSex":"XY","lifeStatus":"ALIVE",
                            "annotationSets":[{"id":"demo2","variableSetId": "individual_private_details",
                                                "annotations": {"full_name": "John Smith", "age": 60,
                                                                "gender": "MALE","hpo": ["HP:0000118", "HP:0000220"]}}]})

    # You can add the annotations during the individual creation step or later through the update method
    individuals.update(query_id="HG00101",
                    data={"sex":"MALE","karyotypicSex":"XY","lifeStatus":"ALIVE","ethnicity":"AFR",
                            "annotationSets":[{"id":"demo2","variableSetId": "individual_private_details",
                                                "annotations": {"full_name": "Bob Smith", "age": 62,
                                                                "gender": "MALE","hpo": ["HP:0000118", "HP:0000220"]}}]})

    # Search Individuals
    in_res = individuals.search(study="1kG_phase3",limit=15)

    for r in in_res.results():
        print(f" {r['id']}, {get_individual_details(r).unwrap_or('NA')}")
    # Note that we have not annotated the other individuals so they have no annotations yet