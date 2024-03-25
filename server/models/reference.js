const mongoose = require('mongoose');

// MongoDB Schema and Model setup
const referenceSchema = new mongoose.Schema({
    pubmedId: String,
    DOI: String,
    year: { type: Number, required: true }
});
  
// Create Mongoose Model
const ReferenceData = mongoose.model('references', referenceSchema);

module.exports = ReferenceData;
  