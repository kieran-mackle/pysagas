## v0.3.0 (2022-12-06)

### Feat

- **__init__.py**: expose key classes to top namespace
- **test_cart3d.py**: code running to completion
- assign FlowState to Cells
- **flow.py**: created GasState and FlowState classes
- **Vector**: added unit vector property
- **utilities.py**: added ParaView macro to process components.i.plt
- **geometry.py**: added class methods for Vector and Cell objects

### Fix

- **flow.py**: circular import
- fixed circular import error
- **utilities.py**: fixed dimensionality for parameters

## v0.2.0 (2022-11-30)

### Feat

- implemented use of Cells
- **Cell**: added sensitivity methods to cell
- **geometry.py**: implemented area and normal methods for Cell
- **geometry.py**: implemented operation methods for Vector
- **utilities.py**: added force sensitivity utilities
- **geometry.py**: added vec property to Vector object
- **geometry.py**: added geometry module
- added some utilities
- **flow.py**: added flow state object

### Refactor

- code tidy up
