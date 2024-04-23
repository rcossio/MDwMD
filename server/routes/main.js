const express = require('express');
const router = express.Router();
const verifyToken = require('../middleware/verifyToken');
const path = require('path');

//Public routes
router.get('/', (req, res) => {
  res.sendFile(path.resolve('public', 'browse.html'));
});

router.get('/browse-proteins', (req, res) => {
  res.sendFile(path.resolve('public', 'browse_proteins.html'));
});

router.get('/login', (req, res) => {
  res.sendFile(path.resolve('public', 'login.html'));
});

router.get('/unauthorized', (req, res) => {
  res.sendFile(path.resolve('public', 'unauthorized.html'));
});

router.get('/privacy', (req, res) => {
  res.sendFile(path.resolve('public', 'privacy.html'));
});

router.get('/contact', (req, res) => {
  res.sendFile(path.resolve('public', 'contact.html'));
});

//Protected routes
router.get('/new_experiment',verifyToken, (req, res) => {
  res.sendFile(path.resolve('public', 'new_experiment.html'));
});

router.get('/new_protein',verifyToken, (req, res) => {
  res.sendFile(path.resolve('public', 'new_protein.html'));
});

router.get('/update/experiment/:id',verifyToken, (req, res) => {
  res.redirect('/update_experiment.html?id='+req.params.id);
});

router.get('/update/protein/:id',verifyToken, (req, res) => {
  res.redirect('/update_protein.html?id='+req.params.id);
});

module.exports = router;
